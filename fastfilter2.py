#!/usr/bin/env python3

"""
FastFilter: Threaded paired-end FASTQ filter for STAR
Author: Lucas Monteiro | PI: Margarida Gama-Carvalho
RNA Systems Biology Lab, BioISI, University of Lisbon

Features:
- Filters sequences by length, quality, homopolymer runs, and N characters
- Supports .fastq and .fastq.gz input
- Produces .fastq.gz output compatible with STAR
- Multi-threaded with stacked progress bars per file
- Dry-run option
- Continuous TSV summary updates
- Logging instead of print statements
"""

import argparse
import itertools
import threading
import gzip
import sys
import time
import csv
import logging
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

from Bio.SeqIO.QualityIO import FastqPhredIterator
from tqdm import tqdm

# --------------------------- Default Parameters --------------------------- #
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30

# --------------------------- Logging Setup --------------------------- #
logging.basicConfig(
    format="%(asctime)s | %(levelname)s | %(message)s",
    level=logging.INFO,
    datefmt="%Y-%m-%d %H:%M:%S"
)

# --------------------------- Sequence Analysis --------------------------- #
def analyze_sequence(record, min_len=MIN_LENGTH_DEFAULT, min_score=MIN_SCORE_DEFAULT, homopolymer_coeff=HOMOPOLYMER_COEFF_DEFAULT):
    """Check if a sequence passes length, quality, homopolymer, and N filters."""
    seq = str(record.seq)
    quals = record.letter_annotations.get("phred_quality", [])

    if len(seq) < min_len:
        return False
    mean_q = sum(quals)/len(quals) if quals else 0
    if mean_q < min_score:
        return False
    max_homopolymer = max((len(list(g)) for _, g in itertools.groupby(seq)), default=0)
    if max_homopolymer > homopolymer_coeff:
        return False
    if "N" in seq or "." in seq:
        return False
    return True

# --------------------------- File Handling --------------------------- #
def open_input(file_path):
    """Open FASTQ or gzipped FASTQ file for reading."""
    return gzip.open(file_path, "rt") if str(file_path).endswith(".gz") else open(file_path, "r")

def open_output(file_path: Path, dry_run=False):
    """Open output FASTQ as gzipped for STAR; /dev/null for dry-run."""
    if dry_run:
        return open("/dev/null", "w")
    return gzip.open(file_path.with_suffix(".fastq.gz"), "wt")

# --------------------------- File Filtering --------------------------- #
def parse_file(r1_file, r2_file, output_dir, min_len, min_score, homopolymer_coeff, dry_run, position, results_list, lock):
    """Filter a single R1/R2 FASTQ pair with stacked progress bar."""
    sample_name = r1_file.stem.replace("_R1","")
    logging.info(f"Processing sample: {sample_name}")

    total_reads = 0
    passed_pairs = 0
    passed_r1 = 0
    passed_r2 = 0

    out_r1 = open_output(output_dir / f"{r1_file.stem}_FILTERED", dry_run)
    out_r2 = open_output(output_dir / f"{r2_file.stem}_FILTERED", dry_run)

    with open_input(r1_file) as f1, open_input(r2_file) as f2, out_r1, out_r2:
        r1_iter = FastqPhredIterator(f1)
        r2_iter = FastqPhredIterator(f2)
        iterator = itertools.zip_longest(r1_iter, r2_iter)
        iterator = tqdm(iterator, desc=sample_name, unit=" reads", position=position, leave=True)

        for rec1, rec2 in iterator:
            if rec1 is None or rec2 is None:
                raise RuntimeError(f"Unequal reads in {r1_file.name} and {r2_file.name}")
            if rec1.id != rec2.id:
                raise RuntimeError(f"Read ID mismatch:\n{rec1.id}\n{rec2.id}")

            total_reads += 1
            pass1 = analyze_sequence(rec1, min_len, min_score, homopolymer_coeff)
            pass2 = analyze_sequence(rec2, min_len, min_score, homopolymer_coeff)
            passed_r1 += pass1
            passed_r2 += pass2

            if pass1 and pass2:
                out_r1.write(rec1.format("fastq"))
                out_r2.write(rec2.format("fastq"))
                passed_pairs += 1

    result = {
        "sample": sample_name,
        "input_reads": total_reads,
        "passed_pairs": passed_pairs,
        "passed_r1": passed_r1,
        "passed_r2": passed_r2
    }

    with lock:
        results_list.append(result)

    logging.info(f"Finished sample: {sample_name} | {passed_pairs}/{total_reads} pairs passed")
    return result

# --------------------------- Summary Handling --------------------------- #
def write_summary(results, output_dir):
    """Write TSV summary of all samples."""
    summary_path = output_dir / "filtering_summary.tsv"
    with open(summary_path, "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["Sample","Input_Reads","Passed_Pairs","Passed_R1","Passed_R2","Percent_Pairs_Passed"])
        for r in results:
            percent = (r["passed_pairs"]/r["input_reads"]*100) if r["input_reads"]>0 else 0
            writer.writerow([r["sample"], r["input_reads"], r["passed_pairs"], r["passed_r1"], r["passed_r2"], f"{percent:.2f}"])
    logging.info(f"Summary written to: {summary_path}")
    return summary_path

# --------------------------- Main --------------------------- #
def main():
    parser = argparse.ArgumentParser(description="FastFilter: threaded paired-end FASTQ filter for STAR")
    parser.add_argument("-l","--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p","--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s","--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i","--sequences-dir", type=str, required=True)
    parser.add_argument("-o","--output-dir", type=str)
    parser.add_argument("-j","--cpus", type=int, default=1)
    parser.add_argument("-d","--dryrun", action="store_true")
    args = parser.parse_args()

    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Match R1/R2 files
    r1_files = sorted(seq_dir.glob("*_R1*.fastq*"))
    r2_files = sorted(seq_dir.glob("*_R2*.fastq*"))
    r1_dict = {f.name.replace("_R1","_"): f for f in r1_files}
    r2_dict = {f.name.replace("_R2","_"): f for f in r2_files}

    missing_r2 = set(r1_dict) - set(r2_dict)
    missing_r1 = set(r2_dict) - set(r1_dict)
    if missing_r2: sys.exit(f"Missing R2 files for: {missing_r2}")
    if missing_r1: sys.exit(f"Missing R1 files for: {missing_r1}")

    paired_files = [(r1_dict[k], r2_dict[k]) for k in sorted(r1_dict.keys())]
    logging.info(f"Detected {len(paired_files)} paired-end samples")

    start_time = time.time()
    results_list = []
    summary_lock = threading.Lock()

    # Threaded execution with stacked progress bars
    with ThreadPoolExecutor(max_workers=args.cpus) as executor:
        futures = []
        for idx, (r1,r2) in enumerate(paired_files):
            futures.append(executor.submit(parse_file, r1, r2, output_dir,
                                           args.minlen, args.min_score, args.homopolymerlen,
                                           args.dryrun, idx, results_list, summary_lock))
        for f in as_completed(futures):
            f.result()  # raise exceptions if any

    write_summary(results_list, output_dir)
    elapsed = (time.time() - start_time)/60
    logging.info(f"Filtering completed in {elapsed:.2f} minutes")

if __name__ == "__main__":
    main()
