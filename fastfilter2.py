#!/usr/bin/env python3.9
"""
FastFilter2: Parallel paired-end FASTQ filter
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
Lab: RNA Systems Biology Lab, BioISI, University of Lisbon
"""

import argparse
import gzip
from pathlib import Path
import csv
import time
from itertools import groupby
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio import SeqIO
from multiprocessing import Pool, cpu_count

# --------------------------- Defaults --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
CHUNK_SIZE = 100_000  # reads per chunk for parallel processing

# --------------------------- Sequence Analysis --------------------------- #
def analyze_sequence(record, min_len, min_score, homopolymer_coeff):
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])
    if len(seq) < min_len:
        return False
    mean_q = sum(quals)/len(quals) if quals else 0
    if mean_q < min_score:
        return False
    if max((len(list(g)) for _, g in groupby(seq)), default=0) > homopolymer_coeff:
        return False
    if "N" in seq or "." in seq:
        return False
    return True

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path, dry_run=False):
    if dry_run:
        return open("/dev/null","w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- Chunking --------------------------- #
def chunked_reads(r1_iter, r2_iter, chunk_size=CHUNK_SIZE):
    chunk = []
    for rec1, rec2 in zip(r1_iter, r2_iter):
        chunk.append((rec1, rec2))
        if len(chunk) >= chunk_size:
            yield chunk
            chunk = []
    if chunk:
        yield chunk

# --------------------------- Processing --------------------------- #
def process_chunk(args):
    chunk, min_len, min_score, homopolymer_coeff = args
    filtered_r1, filtered_r2 = [], []
    passed_pairs = 0
    for r1, r2 in chunk:
        if analyze_sequence(r1, min_len, min_score, homopolymer_coeff) and \
           analyze_sequence(r2, min_len, min_score, homopolymer_coeff):
            filtered_r1.append(r1)
            filtered_r2.append(r2)
            passed_pairs += 1
    return filtered_r1, filtered_r2, len(chunk), passed_pairs

# --------------------------- File Processing --------------------------- #
def process_paired_file(r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, dry_run, threads):
    sample_name = r1_file.stem.replace("_R1","")
    out_r1 = open_output(output_dir / f"{r1_file.stem}_FILTERED", dry_run)
    out_r2 = open_output(output_dir / f"{r2_file.stem}_FILTERED", dry_run)

    total_reads = 0
    passed_pairs = 0

    with open_input(r1_file) as f1, open_input(r2_file) as f2, out_r1, out_r2:
        r1_iter = FastqPhredIterator(f1)
        r2_iter = FastqPhredIterator(f2)
        chunks = list(chunked_reads(r1_iter, r2_iter))

        # Prepare arguments for parallel map
        args_list = [(chunk, min_len, min_score, homopolymer_coeff) for chunk in chunks]

        with Pool(processes=threads) as pool:
            for filtered_r1, filtered_r2, chunk_total, chunk_passed in pool.map(process_chunk, args_list):
                SeqIO.write(filtered_r1, out_r1, "fastq")
                SeqIO.write(filtered_r2, out_r2, "fastq")
                total_reads += chunk_total
                passed_pairs += chunk_passed

    print(f"{sample_name}: Finished | Total reads: {total_reads} | Passed pairs: {passed_pairs}")
    return {"sample": sample_name, "total_reads": total_reads, "passed_pairs": passed_pairs}

# --------------------------- Summary --------------------------- #
def write_summary(results, output_dir):
    summary_path = output_dir / "filtering_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample","Total_Reads","Passed_Pairs","Percent_Passed"])
        for r in results:
            percent = (r["passed_pairs"]/r["total_reads"]*100) if r["total_reads"]>0 else 0
            writer.writerow([r["sample"], r["total_reads"], r["passed_pairs"], f"{percent:.2f}"])
    print(f"Summary written to: {summary_path}")

# --------------------------- Main --------------------------- #
def main():
    parser = argparse.ArgumentParser(description="FastFilter2: parallel paired-end FASTQ filter")
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-t","--threads", type=int, default=cpu_count())
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect paired files
    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}
    missing_r2 = set(r1_dict) - set(r2_dict)
    missing_r1 = set(r2_dict) - set(r1_dict)
    if missing_r2: raise RuntimeError(f"Missing R2 files for: {missing_r2}")
    if missing_r1: raise RuntimeError(f"Missing R1 files for: {missing_r1}")
    paired_files = [(r1_dict[k], r2_dict[k]) for k in sorted(r1_dict.keys())]

    start_time = time.time()
    results = []
    for r1,r2 in paired_files:
        res = process_paired_file(r1, r2, output_dir, args.minlen, args.min_score, args.homopolymerlen, args.dryrun, args.threads)
        results.append(res)

    write_summary(results, output_dir)
    elapsed = (time.time()-start_time)/60
    print(f"Filtering completed in {elapsed:.2f} minutes")

if __name__ == "__main__":
    main()
