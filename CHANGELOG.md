# FastFilter2 Upgrade Guide

This document outlines the changes, improvements, and migration notes for users moving from **FastFilter v1** to **FastFilter2**.

---

## **Version Overview**

| Feature | FastFilter v1 | FastFilter2 |
|---------|---------------|-------------|
| Input file support | `.fastq` | `.fastq` & `.fastq.gz` |
| Processing | Single-threaded or multiprocessing | Multi-threaded via `ThreadPoolExecutor` |
| Paired-end support | Yes, with R1/R2 parsing | Yes, with strict ID matching & progress bars |
| Output | FASTQ files | gzipped FASTQ (`.fastq.gz`) compatible with STAR |
| Logging | Print statements | Thread-safe logging with timestamps |
| Summary | CSV & TXT files | Continuous TSV summary with detailed stats |
| Dry-run | Yes | Yes, outputs skipped using `/dev/null` |
| Progress monitoring | None | Stacked progress bars per sample with `tqdm` |

---

## **Key Improvements in FastFilter2**

1. **Multi-threading:**  
   - Each sample pair is processed in a separate thread, significantly reducing runtime on multi-core systems.  
   - Progress bars are stacked per sample for visual monitoring.

2. **Gzipped Input/Output:**  
   - Supports `.fastq.gz` files directly.  
   - Outputs are gzipped for compatibility with STAR and storage efficiency.

3. **Enhanced Logging:**  
   - Uses `logging` instead of `print` for professional timestamped messages.  
   - Errors, mismatched read IDs, or missing R1/R2 files are clearly reported.

4. **Strict Paired-End Validation:**  
   - Checks that R1 and R2 have equal read counts and matching IDs.  
   - Raises errors immediately if mismatches are found.

5. **Continuous TSV Summary:**  
   - Generates `filtering_summary.tsv` during processing.  
   - Includes total reads, passed pairs, reads passing individually (R1/R2), and percent passed.

6. **Dry-Run Mode:**  
   - Allows testing without writing any output files, using `/dev/null` as a sink.

7. **Simplified Argument Handling:**  
   - Optional arguments for minimum length, quality score, homopolymer length, and CPU threads.  
   - Automatic creation of output directories.

---

## **Migration Notes**

- **Output File Names:**  
  - Old: `<sample>_FILTERED.fastq`  
  - New: `<sample>_R1_FILTERED.fastq.gz` and `<sample>_R2_FILTERED.fastq.gz`

- **Summary Format:**  
  - Old: CSV with separate R1/R2 exclusion reasons  
  - New: TSV with passed pairs and percent passed per sample

- **Threading & Progress Bars:**  
  - Previous scripts could hang if multiprocessing failed; FastFilter2 uses threads with robust exception handling.

- **Dry-Run Behavior:**  
  - Outputs are sent to `/dev/null` instead of skipping entirely.

- **Error Handling:**  
  - Mismatched R1/R2 files or read counts will terminate processing with clear logs, instead of silently skipping.

---

## **Example Migration Command**

**Old command (v1):**

python fastfilter.py -i /data/project/cutadapt -o /data/project/fastfilter -l 30  

**New command (v2):**

python FastFilter2.py -i /data/project/cutadapt -o /data/project/fastfilter -l 30 -s 30 -p 25 -j 4  

- Adds explicit threading (`-j 4`)  
- Adds Phred score filter (`-s 30`)  
- Adds homopolymer max length (`-p 25`)  

---

## **Recommendations**

1. Test FastFilter2 on a small dataset using `--dryrun` before full-scale filtering.  
2. Validate STAR compatibility with gzipped outputs.  
3. Review `filtering_summary.tsv` to ensure thresholds match your experimental design.  

---

## **References**

- FastFilter original implementation: [fastfilter README](https://github.com/GamaPintoLab/fastfilter/blob/main/README.md)
- STAR aligner: https://github.com/alexdobin/STAR  

---
