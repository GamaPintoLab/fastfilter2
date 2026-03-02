#!/usr/bin/env python3.9

"""
fastfilter_pe.py - Paired-end FASTQ filtering for STAR
FIXED VERSION:
- No file overwriting
- Proper batch writing
- Writes directly to gzip
- Safe multiprocessing
"""

import argparse
import gzip
import multiprocessing
import sys
from pathlib import Path
from Bio.SeqIO.QualityIO import FastqPhredIterator
from Bio import SeqIO
from tqdm import tqdm

# Defaults
MIN_LENGTH_DEFAULT = 25
HOMOPOLYMER_COEFF_DEFAULT = 25
MIN_SCORE_DEFAULT = 30
WRITE_BATCH_SIZE = 1000

# Globals
min_seq_len = MIN_LENGTH_DEFAULT
homopolymer_coeff = HOMOPOLYMER_COEFF_DEFAULT
min_score = MIN_SCORE_DEFAULT
num_cpus = 1
seq_dir = None
output_dir = None
dryrun = False


# ---------------------
# Argument parsing
# ---------------------
def parse_arguments():
    global min_seq_len, homopolymer_coeff, min_score
    global seq_dir, output_dir, num_cpus, dryrun

    parser = argparse.ArgumentParser(description="Paired-end FASTQ filtering for STAR.")
    parser.add_argument("-l", "--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p", "--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s", "--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i", "--seq-dir", type=str, required=True)
    parser.add_argument("-o", "--output-dir", type=str)
    parser.add_argument("-j", "--cpus", type=int, default=1)
    parser.add_argument("--dryrun", action="store_true")

    args = parser.parse_args()

    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.min_score
    seq_dir = Path(args.seq_dir)
    output_dir = Path(args.output_dir) if args.output_dir else seq_dir.parent / "fastfilter"
    num_cpus = args.cpus
    dryrun = args.dryrun


# ---------------------
# Homopolymer detection
# ---------------------
def find_homopolymers(seq):
    # Faster and correct contiguous detection
    for nt in "ATGC":
        if nt * homopolymer_coeff in seq:
            return True
    return False


# ---------------------
# Sequence filter
# ---------------------
def filter_sequence(record):
    seq = str(record.seq)
    qual = record.letter_annotations.get("phred_quality", [])

    if len(seq) < min_seq_len:
        return False
    if "N" in seq or "." in seq:
        return False
    if find_homopolymers(seq):
        return False
    if not qual:
        return False
    if sum(qual) / len(qual) < min_score:
        return False

    return True


# ---------------------
# Process paired-end file
# ---------------------
def process_pair(r1_path: Path, r2_path: Path, position: int):
    pair_name = r1_path.stem.replace("_R1", "")

    r1_out_path = output_dir / f"{pair_name}_R1_FILTERED.fastq.gz"
    r2_out_path = output_dir / f"{pair_name}_R2_FILTERED.fastq.gz"

    total_reads = 0
    good_reads = 0

    r1_buffer = []
    r2_buffer = []

    # Open output once (no overwriting issue anymore)
    if not dryrun:
        r1_out = gzip.open(r1_out_path, "wt")
        r2_out = gzip.open(r2_out_path, "wt")
    else:
        r1_out = None
        r2_out = None

    with open(r1_path, "r") as r1_in, open(r2_path, "r") as r2_in:
        r1_iterator = FastqPhredIterator(r1_in)
        r2_iterator = FastqPhredIterator(r2_in)

        for rec1, rec2 in tqdm(
            zip(r1_iterator, r2_iterator),
            desc=pair_name,
            unit="reads",
            position=position,
            leave=True,
            dynamic_ncols=True
        ):
            total_reads += 1

            keep1 = filter_sequence(rec1)
            keep2 = filter_sequence(rec2)

            if keep1 and keep2:
                good_reads += 1

                if not dryrun:
                    r1_buffer.append(rec1)
                    r2_buffer.append(rec2)

                # Batch write
                if not dryrun and len(r1_buffer) >= WRITE_BATCH_SIZE:
                    SeqIO.write(r1_buffer, r1_out, "fastq")
                    SeqIO.write(r2_buffer, r2_out, "fastq")
                    r1_buffer.clear()
                    r2_buffer.clear()

    # Write remaining
    if not dryrun and r1_buffer:
        SeqIO.write(r1_buffer, r1_out, "fastq")
        SeqIO.write(r2_buffer, r2_out, "fastq")

    if not dryrun:
        r1_out.close()
        r2_out.close()

    return {
        "file": pair_name,
        "total_reads": total_reads,
        "good_reads": good_reads
    }


# ---------------------
# Main workflow
# ---------------------
def main():
    parse_arguments()

    if not seq_dir.is_dir():
        print(f"Input directory not found: {seq_dir}")
        sys.exit(1)

    output_dir.mkdir(parents=True, exist_ok=True)

    r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))

    if not r1_files or len(r1_files) != len(r2_files):
        print("Error: R1/R2 mismatch.")
        sys.exit(1)

    pairs = list(zip(r1_files, r2_files))
    positions = list(range(len(pairs)))

    pool_args = [(r1, r2, pos) for (r1, r2), pos in zip(pairs, positions)]

    with multiprocessing.Pool(processes=num_cpus) as pool:
        results = pool.starmap(process_pair, pool_args)

    # Summary
    summary_file = output_dir / "fastfilter_summary.csv"

    with open(summary_file, "w") as f:
        f.write("file,total_reads,good_reads,pass_rate_pct\n")
        for r in results:
            rate = (r["good_reads"] / r["total_reads"] * 100) if r["total_reads"] else 0
            f.write(f"{r['file']},{r['total_reads']},{r['good_reads']},{rate:.2f}\n")

    print(f"\nFiltering complete.")
    print(f"Summary written to: {summary_file}")


if __name__ == "__main__":
    main()
