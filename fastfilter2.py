#!/usr/bin/env python3.9

"""
FastFilter2: High-performance threaded FASTQ filter
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon

Features:
- Filters sequences by length, quality, homopolymer runs, and N characters
- Supports .fastq and .fastq.gz input
- Produces .fastq.gz output compatible with STAR
- Processes one paired file at a time with threads for speed
- Dry-run option
- TSV summary of all samples
- Stacked progress bar per sample
"""

import argparse
import gzip
import itertools
import csv
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
from Bio.SeqIO.QualityIO import FastqPhredIterator
from tqdm import tqdm

# --------------------------- Defaults --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
CHUNK_SIZE = 100_000  # batch size for threads

# --------------------------- Analysis --------------------------- #
def analyze_sequence(record, min_len, min_score, homopolymer_coeff):
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])
    if len(seq) < min_len:
        return False
    mean_q = sum(quals)/len(quals) if quals else 0
    if mean_q < min_score:
        return False
    if max((len(list(g)) for _, g in itertools.groupby(seq)), default=0) > homopolymer_coeff:
        return False
    if "N" in seq or "." in seq:
        return False
    return True

def process_chunk(chunk, min_len, min_score, homopolymer_coeff):
    filtered_r1 = []
    filtered_r2 = []
    stats = {"passed_pairs":0, "passed_r1":0, "passed_r2":0, "total":0}
    for r1, r2 in chunk:
        stats["total"] += 1
        p1 = analyze_sequence(r1, min_len, min_score, homopolymer_coeff)
        p2 = analyze_sequence(r2, min_len, min_score, homopolymer_coeff)
        stats["passed_r1"] += p1
        stats["passed_r2"] += p2
        if p1 and p2:
            stats["passed_pairs"] += 1
            filtered_r1.append(r1.format("fastq"))
            filtered_r2.append(r2.format("fastq"))
    return filtered_r1, filtered_r2, stats

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path, dry_run=False):
    if dry_run:
        return open("/dev/null","w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- Single File --------------------------- #
def process_pair_file(r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, dry_run, threads):
    sample = r1_file.stem.replace("_R1","")
    out_r1 = open_output(output_dir / f"{r1_file.stem}_FILTERED", dry_run)
    out_r2 = open_output(output_dir / f"{r2_file.stem}_FILTERED", dry_run)
    
    total_reads = 0
    passed_pairs = 0
    passed_r1 = 0
    passed_r2 = 0

    with open_input(r1_file) as f1, open_input(r2_file) as f2, out_r1, out_r2:
        r1_iter = FastqPhredIterator(f1)
        r2_iter = FastqPhredIterator(f2)
        batch = []
        futures = []
        executor = ThreadPoolExecutor(max_workers=threads) if threads > 1 else None
        tqdm_bar = tqdm(unit=" reads", desc=sample)

        for rec1, rec2 in itertools.zip_longest(r1_iter, r2_iter):
            if rec1 is None or rec2 is None:
                raise RuntimeError(f"Unequal reads in {r1_file.name} and {r2_file.name}")
            if rec1.id != rec2.id:
                raise RuntimeError(f"Read ID mismatch:\n{rec1.id}\n{rec2.id}")

            batch.append((rec1, rec2))
            if len(batch) >= CHUNK_SIZE:
                if executor:
                    futures.append(executor.submit(process_chunk, batch, min_len, min_score, homopolymer_coeff))
                else:
                    fr1, fr2, stats = process_chunk(batch, min_len, min_score, homopolymer_coeff)
                    out_r1.writelines(fr1)
                    out_r2.writelines(fr2)
                    total_reads += stats["total"]
                    passed_pairs += stats["passed_pairs"]
                    passed_r1 += stats["passed_r1"]
                    passed_r2 += stats["passed_r2"]
                    tqdm_bar.update(stats["total"])
                batch = []

        if batch:
            if executor:
                futures.append(executor.submit(process_chunk, batch, min_len, min_score, homopolymer_coeff))
            else:
                fr1, fr2, stats = process_chunk(batch, min_len, min_score, homopolymer_coeff)
                out_r1.writelines(fr1)
                out_r2.writelines(fr2)
                total_reads += stats["total"]
                passed_pairs += stats["passed_pairs"]
                passed_r1 += stats["passed_r1"]
                passed_r2 += stats["passed_r2"]
                tqdm_bar.update(stats["total"])

        if executor:
            for fut in as_completed(futures):
                fr1, fr2, stats = fut.result()
                out_r1.writelines(fr1)
                out_r2.writelines(fr2)
                total_reads += stats["total"]
                passed_pairs += stats["passed_pairs"]
                passed_r1 += stats["passed_r1"]
                passed_r2 += stats["passed_r2"]
                tqdm_bar.update(stats["total"])
            executor.shutdown(wait=True)

        tqdm_bar.close()

    print(f"Sample {sample}: Total reads={total_reads}, Passed pairs={passed_pairs}, "
          f"R1 passed={passed_r1}, R2 passed={passed_r2}")
    
    return {
        "sample": sample,
        "input_reads": total_reads,
        "passed_pairs": passed_pairs,
        "passed_r1": passed_r1,
        "passed_r2": passed_r2
    }

# --------------------------- Summary --------------------------- #
def write_summary(results, output_dir):
    path = output_dir / "filtering_summary.tsv"
    with open(path,"w",newline="") as f:
        writer = csv.writer(f,delimiter="\t")
        writer.writerow(["Sample","Input_Reads","Passed_Pairs","Passed_R1","Passed_R2","Percent_Pairs_Passed"])
        for r in results:
            pct = r["passed_pairs"]/r["input_reads"]*100 if r["input_reads"]>0 else 0
            writer.writerow([r["sample"],r["input_reads"],r["passed_pairs"],r["passed_r1"],r["passed_r2"],f"{pct:.2f}"])
    print(f"Summary written to: {path}")

# --------------------------- Main --------------------------- #
def main():
    parser = argparse.ArgumentParser(description="FastFilter2 threaded FASTQ filter")
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-j","--threads", type=int, default=1)
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Match R1/R2
    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}
    missing_r2 = set(r1_dict) - set(r2_dict)
    missing_r1 = set(r2_dict) - set(r1_dict)
    if missing_r2: raise RuntimeError(f"Missing R2 for: {missing_r2}")
    if missing_r1: raise RuntimeError(f"Missing R1 for: {missing_r1}")

    paired_files = [(r1_dict[k], r2_dict[k]) for k in sorted(r1_dict.keys())]

    results = []
    for r1,r2 in paired_files:
        res = process_pair_file(r1,r2,output_dir,args.minlen,args.min_score,args.homopolymerlen,args.dryrun,args.threads)
        results.append(res)

    write_summary(results, output_dir)

if __name__ == "__main__":
    main()
