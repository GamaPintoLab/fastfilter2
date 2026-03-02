#!/usr/bin/env python3.9
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
    parser.add_argument("-l", "--minlen", type=int, default=MIN_LENGTH_DEFAULT)
    parser.add_argument("-p", "--homopolymerlen", type=int, default=HOMOPOLYMER_COEFF_DEFAULT)
    parser.add_argument("-s", "--min-score", type=int, default=MIN_SCORE_DEFAULT)
    parser.add_argument("-i", "--seq-dir", type=str, required=True)
    parser.add_argument("-o", "--output-dir", type=str)
    parser.add_argument("-j", "--cpus", type=int, default=1)
    args = parser.parse_args()

    min_seq_len = args.minlen
    homopolymer_coeff = args.homopolymerlen
    min_score = args.min_score
    seq_dir = Path(args.seq_dir)
    output_dir = Path(args.output_dir) if args.output_dir else None
    num_cpus = args.cpus

# ---------------------
# Homopolymer detection
# ---------------------
def has_homopolymer(seq):
    return any(nt * homopolymer_coeff in seq for nt in "ATGC")

# ---------------------
# Sequence filter
# ---------------------
def filter_sequence(record):
    seq = record.seq
    qual = record.letter_annotations.get("phred_quality", [])
    avg_qual = sum(qual)/len(qual) if qual else 0
    if len(seq) < min_seq_len: return False
    if "N" in seq or "." in seq: return False
    if has_homopolymer(seq): return False
    if avg_qual < min_score: return False
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
            r1_iter = FastqPhredIterator(r1_in)
            r2_iter = FastqPhredIterator(r2_in)

            for rec1, rec2 in tqdm(zip(r1_iter, r2_iter),
                                   desc=f"{pair_name}",
                                   unit="reads",
                                   ncols=80,
                                   leave=True,
                                   unit_scale=True):
                total_reads += 1
                if filter_sequence(rec1) and filter_sequence(rec2):
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

    if not seq_dir.is_dir():
        print(f"Input sequence directory not found: {seq_dir}")
        sys.exit(1)

    output_dir = output_dir or seq_dir.parent / "fastfilter"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Paired-end files
    r1_files = sorted(seq_dir.glob("*_R1_*.fastq"))
    r2_files = sorted(seq_dir.glob("*_R2_*.fastq"))

    if not r1_files or not r2_files or len(r1_files) != len(r2_files):
        print("Error: Paired-end R1/R2 files missing or mismatch.")
        sys.exit(1)

    pairs = list(zip(r1_files, r2_files))

    # Multiprocessing
    with multiprocessing.Pool(processes=num_cpus) as pool:
        results = pool.starmap(process_pair, pairs)

    # Summary CSV
    summary_file = output_dir / "fastfilter_summary.csv"
    with open(summary_file, "w") as f:
        f.write("file,total_reads,good_reads,pass_rate_pct\n")
        for r in results:
            rate = (r["good_reads"] / r["total_reads"] * 100) if r["total_reads"] else 0
            f.write(f"{r['file']},{r['total_reads']},{r['good_reads']},{rate:.1f}\n")

    print(f"Filtering complete. Summary saved to {summary_file}")

if __name__ == "__main__":
    main()
