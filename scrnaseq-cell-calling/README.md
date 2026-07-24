# scRNA-seq cell calling

A synthetic, multiplexed scRNA-seq QC demo for cell calling, ambient-RNA
estimation, and per-sample summaries. The count matrix uses generic probe and
sample barcodes.

## Run

```bash
pip install -r requirements.txt
make scrna-cell-calling
```

## What it does

- Simulates probe-level UMI counts for four multiplexed samples.
- Adds ambient background and empty barcodes.
- Calls cells from the barcode-rank knee.
- Scores calls against planted truth.
- Reports cell count, median UMIs, and detected probes per sample.

The knee is the point of maximum curvature on a smoothed log-log barcode-rank
curve. Barcodes at or above the corresponding UMI count are called as cells.

![scRNA-seq cell-calling QC](assets/scrnaseq_qc.png)

## Files

```text
generate_data.py   synthetic multiplexed scRNA-seq counts and empty barcodes
analyze.py         cell calling, ambient profile, and per-sample QC
plots.py           barcode rank, UMI distributions, and ambient profile
```
