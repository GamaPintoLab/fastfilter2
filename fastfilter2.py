#!/usr/bin/env python3.9

"""
fastfilter_pe.py - Paired-end FASTQ filtering for STAR
Author: Updated for paired-end filtering

Filters paired-end FASTQ files:
- Removes sequences that are too short, contain Ns, dots, or long homopolymers.
- Applies min average quality filter.
- Writes filtered R1 and R2 reads to FASTQ.GZ outputs for STAR alignment.
- Generates a summary CSV with total vs good reads.
"""

import argparse
import gzip
import multiprocessing
import sys
import time
from pathlib import Path

from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio import SeqIO
from tqdm import tqdm

# Defaults
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
PROJECTS_DIR_DEFAULT = Path("/data/working_directory/projects")

# Globals
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
num_cpus = 1
seq_dir = None
output_dir = None

# ---------------------
# Argument parsing
# ---------------------
def parse_arguments():
    global min_seq_len, homopolymer_coeff, min_score, seq_dir, output_dir, num_cpus

    parser = argparse.ArgumentParser(description="Paired-end FASTQ filtering for STAR.")
    parser.add_argument("-l", "--minlen", type=int, default=MIN_LENGTH_DEFAULT,
                        help="Minimum sequence length")
    parser.add_argument("-p", "--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT,
                        help="Length threshold for homopolymers")
    parser.add_argument("-s", "--min-score", type=int, default=MIN_SCORE_DEFAULT,
                        help="Minimum average Phred quality score")
    parser.add_argument("-i", "--seq-dir", type=str, help="Input directory containing cutadapt FASTQ files")
    parser.add_argument("-o", "--output-dir", type=str, help="Directory to save filtered FASTQ.gz files")
    parser.add_argument("-j", "--cpus", type=int, default=1, help="Number of parallel CPUs")

    args = parser.parse_args()
    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.min_score
    seq_dir = Path(args.seq_dir) if args.seq_dir else None
    output_dir = Path(args.output_dir) if args.output_dir else None
    num_cpus = args.cpus

# ---------------------
# Homopolymer detection
# ---------------------
def find_homopolymers(seq):
    return any(seq.count(nt * homopolymer_coeff) > 0 for nt in "ATGC")

# ---------------------
# Sequence filter
# ---------------------
def filter_sequence(record):
    seq = str(record.seq)
    qual = record.letter_annotations.get("phred_quality", [])
    avg_qual = sum(qual) / len(qual) if qual else 0

    if len(seq) < min_seq_len:
        return False
    if "N" in seq or "." in seq:
        return False
    if find_homopolymers(seq):
        return False
    if avg_qual < min_score:
        return False
    return True

# ---------------------
# Process paired-end file
# ---------------------
def process_pair(r1_path: Path, r2_path: Path):
    pair_name = r1_path.stem.replace("_R1", "")
    r1_out_path = output_dir / f"{pair_name}_R1_FILTERED.fastq.gz"
    r2_out_path = output_dir / f"{pair_name}_R2_FILTERED.fastq.gz"

    total_reads = 0
    good_reads = 0

    with gzip.open(r1_out_path, "wt") as r1_out, gzip.open(r2_out_path, "wt") as r2_out:
        with open(r1_path, "r") as r1_in, open(r2_path, "r") as r2_in:
            r1_iterator = FastqPhredIterator(r1_in)
            r2_iterator = FastqPhredIterator(r2_in)

            for rec1, rec2 in tqdm(zip(r1_iterator, r2_iterator),
                                   desc=f"Filtering {pair_name}", unit="reads"):
                total_reads += 1
                keep1 = filter_sequence(rec1)
                keep2 = filter_sequence(rec2)
                if keep1 and keep2:
                    SeqIO.write(rec1, r1_out, "fastq")
                    SeqIO.write(rec2, r2_out, "fastq")
                    good_reads += 1

    return {"file": pair_name, "total_reads": total_reads, "good_reads": good_reads}

# ---------------------
# Main workflow
# ---------------------
def main():
    global seq_dir, output_dir
    parse_arguments()

    if seq_dir is None or not seq_dir.is_dir():
        print(f"Input sequence directory not found: {seq_dir}")
        sys.exit(1)
    output_dir = output_dir or seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find paired-end FASTQ files
    r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))

    if not r1_files or not r2_files or len(r1_files) != len(r2_files):
        print("Error: Paired-end R1/R2 files not found or mismatch in counts.")
        sys.exit(1)

    pairs = list(zip(r1_files, r2_files))

    # Multiprocessing
    with multiprocessing.Pool(processes=num_cpus) as pool:
        results = pool.starmap(process_pair, pairs)

    # Write summary CSV
    summary_file = output_dir / "fastfilter_summary.csv"
    with open(summary_file, "w") as f:
        f.write("file,total_reads,good_reads,pass_rate_pct\n")
        for r in results:
            rate = (r["good_reads"] / r["total_reads"] * 100) if r["total_reads"] > 0 else 0
            f.write(f"{r['file']},{r['total_reads']},{r['good_reads']},{rate:.1f}\n")

    print(f"Filtering complete. Summary written to {summary_file}")

if __name__ == "__main__":
    main()
