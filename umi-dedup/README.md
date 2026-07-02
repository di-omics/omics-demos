# umi-dedup

Why UMIs matter, in one runnable example. Simulates a PCR-amplified, UMI-tagged
library, then recovers true molecule counts from duplicated, error-containing reads.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make umi
```

## What it does

- Generates ~3,000 true molecules, tags each with a random UMI, PCR-amplifies to reads.
- A small per-base error rate creates near-duplicate UMIs (the realistic hard case).
- Counts three ways: raw reads, naive-unique UMIs, and a **1-mismatch directional collapse**.
- Reports duplication rate and a library **saturation curve**.

## Example output

```
Total reads:           18000
Naive-unique molecules:3528
Collapsed molecules:   2999   (1-mismatch merge)
Duplication rate:      83.3%
```

Naive counting overcounts (3,528) because sequencing errors mint spurious UMIs;
the 1-mismatch collapse recovers 2,999 - within one of the 3,000 truth.

![UMI QC](assets/umi_qc.png)

## Files

```
generate_data.py   simulate UMI-tagged reads with PCR duplication + UMI errors
analyze.py         raw vs naive vs 1-mismatch-collapsed counts, saturation
plots.py           raw-vs-deduped bars + saturation curve
```
