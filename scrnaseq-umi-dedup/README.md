# scRNA-seq UMI deduplication

A runnable example of UMI-aware molecule counting for scRNA-seq. It simulates
PCR-duplicated, UMI-tagged reads and compares naive, directional, and adjacency
collapse against planted molecule truth.

All data is synthetic.

## Run

```bash
pip install -r requirements.txt
make scrna-umi
```

## What it does

- Generates approximately 3,000 captured molecules across synthetic genes.
- Assigns each molecule a random UMI.
- Simulates PCR duplication and UMI sequencing errors.
- Compares naive unique-UMI counting with one-mismatch collapse methods.
- Produces a library-complexity rarefaction curve.

Directional collapse lets a high-count UMI absorb a one-mismatch, lower-count
neighbor when the count ratio supports an error relationship.

![scRNA-seq UMI QC](assets/scrnaseq_umi_qc.png)

## Files

```text
generate_data.py   synthetic scRNA-seq reads with UMI duplication and errors
analyze.py         molecule counting, collapse methods, and saturation
plots.py           per-gene comparison and rarefaction
```
