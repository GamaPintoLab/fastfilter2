#!/usr/bin/env python3.9
"""
fastfilter_paired.py - version 2.0
Author: Gil Poiares-Oliveira <gpo@ciencias.ulisboa.pt>
PI: Margarida Gama-Carvalho <mhcarvalho@ciencias.ulisboa.pt>
RNA Systems Biology Lab, BioISI
(C) 2026
"""

import argparse
import csv
import gzip
import multiprocessing
import sys
import time
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqPhredIterator
from tqdm import tqdm  # progress bar

# Defaults
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
PROJECTS_DIR_DEFAULT = Path("/data/working_directory/projects")

# Globals
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
threads = 1
seq_dir = None
output_dir = None


def parse_arguments():
    global min_seq_len, homopolymer_coeff, min_score, seq_dir, output_dir, threads

    parser = argparse.ArgumentParser(description="Filter paired-end FASTQ files.")
    parser.add_argument("-i", "--sequences-dir", type=str, required=True,
                        help="Directory containing paired FASTQ files (_R1_*.fastq and _R2_*.fastq)")
    parser.add_argument("-o", "--output-dir", type=str,
                        help="Directory for output filtered FASTQ files")
    parser.add_argument("-l", "--minlen", type=int, default=MIN_LENGTH_DEFAULT,
                        help="Minimum sequence length")
    parser.add_argument("-p", "--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT,
                        help="Maximum homopolymer length allowed")
    parser.add_argument("-s", "--min-score", type=int, default=MIN_SCORE_DEFAULT,
                        help="Minimum average quality score")
    parser.add_argument("-j", "--threads", type=int, default=1,
                        help="Number of parallel threads")
    args = parser.parse_args()

    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.min_score
    threads = args.threads
    seq_dir = Path(args.sequences_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)


def find_homopolymers(seq: str) -> dict:
    """Return a dict counting homopolymers of each base."""
    return {
        homopolymer_coeff * "A": int(seq.count("A" * homopolymer_coeff) > 0),
        homopolymer_coeff * "T": int(seq.count("T" * homopolymer_coeff) > 0),
        homopolymer_coeff * "G": int(seq.count("G" * homopolymer_coeff) > 0),
        homopolymer_coeff * "C": int(seq.count("C" * homopolymer_coeff) > 0),
    }


def analyze_sequence(record):
    qual_values = record.letter_annotations["phred_quality"]
    qual_score = sum(qual_values) / len(qual_values) if qual_values else 0
    seq_len = len(record.seq)
    n_count = record.seq.count("N")
    dot_count = record.seq.count(".")
    homopolymers = find_homopolymers(record.seq)
    homopolymers_exist = any(homopolymers.values())

    meets_criteria = (
        seq_len >= min_seq_len and
        n_count == 0 and
        dot_count == 0 and
        not homopolymers_exist and
        qual_score >= min_score
    )

    return {
        "meets_criteria": meets_criteria,
        "length": seq_len,
        "qual_score": qual_score,
        "n_count": n_count,
        "dot_count": dot_count,
        "homopolymers": homopolymers,
        "too_short": seq_len < min_seq_len,
        "found_homopolymer": homopolymers_exist,
        "low_score": qual_score < min_score
    }


def process_pair(r1_file: Path, r2_file: Path) -> dict:
    """Process a single pair of FASTQ files, return summary stats."""
    r1_out = output_dir / (r1_file.stem.replace("_R1", "") + "_FILTERED_R1.fastq")
    r2_out = output_dir / (r2_file.stem.replace("_R2", "") + "_FILTERED_R2.fastq")

    stats = {
        "file": r1_file.stem.replace("_R1", ""),
        "total": 0,
        "good": 0,
        "R1_homopolymers": {},
        "R2_homopolymers": {}
    }

    with open(r1_file, "r") as r1_handle, open(r2_file, "r") as r2_handle, \
         open(r1_out, "w") as r1_filt, open(r2_out, "w") as r2_filt:

        r1_iter = FastqPhredIterator(r1_handle)
        r2_iter = FastqPhredIterator(r2_handle)

        for rec1, rec2 in tqdm(zip(r1_iter, r2_iter), desc=f"Processing {r1_file.name}", unit="reads"):
            stats["total"] += 1
            seq1 = analyze_sequence(rec1)
            seq2 = analyze_sequence(rec2)

            # Update homopolymer counts
            for k in seq1["homopolymers"]:
                stats["R1_homopolymers"][k] = stats["R1_homopolymers"].get(k, 0) + seq1["homopolymers"][k]
                stats["R2_homopolymers"][k] = stats["R2_homopolymers"].get(k, 0) + seq2["homopolymers"][k]

            if seq1["meets_criteria"] and seq2["meets_criteria"]:
                SeqIO.write(rec1, r1_filt, "fastq")
                SeqIO.write(rec2, r2_filt, "fastq")
                stats["good"] += 1

    # Compress the filtered files to .fastq.gz
    for out_file in [r1_out, r2_out]:
        with open(out_file, 'rb') as f_in, gzip.open(str(out_file) + ".gz", 'wb') as f_out:
            f_out.writelines(f_in)
        out_file.unlink()  # remove original uncompressed file

    return stats


def main():
    parse_arguments()

    r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))

    if len(r1_files) != len(r2_files):
        print("Mismatched R1/R2 files!")
        sys.exit(1)

    args_pairs = list(zip(r1_files, r2_files))
    results = []

    # Multiprocessing pool
    with multiprocessing.Pool(threads) as pool:
        for result in pool.starmap(process_pair, args_pairs):
            results.append(result)

    # Save summary CSV
    summary_csv = output_dir / "fastfilter_summary.csv"
    with open(summary_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["file", "total_reads", "good_reads"])
        for r in results:
            writer.writerow([r["file"], r["total"], r["good"]])

    print(f"Filtering completed. Summary saved to {summary_csv}")


if __name__ == "__main__":
    main()
