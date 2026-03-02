#!/usr/bin/env python3.9

"""
fastfilter_pe.py - Paired-end FASTQ filtering for STAR
Author: Updated for speed with final compression

Changes:
- Writes plain FASTQ during filtering (faster than streaming gzip)
- Compresses outputs after filtering
- Fixed tqdm progress bars to stack properly in multiprocessing
"""

import argparse
import gzip
import multiprocessing
import sys
from pathlib import Path
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio import SeqIO
from tqdm import tqdm
import subprocess

# Defaults
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30

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
    parser.add_argument("-o", "--output-dir", type=str, help="Directory to save filtered FASTQ files")
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
def process_pair(r1_path: Path, r2_path: Path, position: int):
    pair_name = r1_path.stem.replace("_R1", "")
    r1_tmp = output_dir / f"{pair_name}_R1_FILTERED.fastq"
    r2_tmp = output_dir / f"{pair_name}_R2_FILTERED.fastq"

    good_reads = 0

    with open(r1_path, "r") as r1_in, open(r2_path, "r") as r2_in, \
         open(r1_tmp, "w") as r1_out, open(r2_tmp, "w") as r2_out:

        r1_iterator = FastqPhredIterator(r1_in)
        r2_iterator = FastqPhredIterator(r2_in)

        # Use tqdm with position for stacked bars in multiprocessing
        for rec1, rec2 in tqdm(zip(r1_iterator, r2_iterator),
                               desc=f"{pair_name}",
                               unit="reads",
                               position=position,
                               leave=True,
                               dynamic_ncols=True):
            keep1 = filter_sequence(rec1)
            keep2 = filter_sequence(rec2)
            if keep1 and keep2:
                SeqIO.write(rec1, r1_out, "fastq")
                SeqIO.write(rec2, r2_out, "fastq")
                good_reads += 1

    # Compress files at the end (much faster)
    for tmp in [r1_tmp, r2_tmp]:
        subprocess.run(["gzip", "-f", str(tmp)])

    return {"file": pair_name, "good_reads": good_reads}

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
    positions = list(range(len(pairs)))  # for stacked progress bars

    # Multiprocessing
    pool_args = [(r1, r2, pos) for (r1, r2), pos in zip(pairs, positions)]
    with multiprocessing.Pool(processes=num_cpus) as pool:
        results = pool.starmap(process_pair, pool_args)

    # Write summary CSV
    summary_file = output_dir / "fastfilter_summary.csv"
    with open(summary_file, "w") as f:
        f.write("file,good_reads\n")
        for r in results:
            f.write(f"{r['file']},{r['good_reads']}\n")

    print(f"Filtering complete. Summary written to {summary_file}")

if __name__ == "__main__":
    main()
