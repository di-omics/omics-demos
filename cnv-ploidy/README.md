# cnv-ploidy

Copy-number variation detection from binned whole-genome coverage, in one
runnable example. Normalizes coverage, corrects for GC bias, and calls gains
and losses with a simple segmentation approach.

All data is synthetic -- no real sequencing data.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make cnv
```

## What it does

- Generates ~1,000 genomic bins across 5 chromosomes with ~100x Poisson coverage.
- Plants three CNV events: a 1.5x gain, a 0.5x loss, and a 2.0x gain.
- Adds mild GC bias (~10% sinusoidal modulation).
- Normalizes by median coverage, corrects GC bias, computes log2 ratio.
- Smooths with a rolling median and calls gains/losses by threshold.
- Compares called segments to planted truth.

## Example output

```
Total bins:       1000
Median coverage:  100.0
Segments called:  3

  chr2:50-80  gain (ratio 1.5) -> RECOVERED
  chr3:100-120 loss (ratio 0.5) -> RECOVERED
  chr5:30-50  gain (ratio 2.0) -> RECOVERED
```

Planted CNVs recovered; note resolution limits at small events.

![CNV QC](assets/cnv_qc.png)

## How it works

GC-bias correction and log2-ratio segmentation (from `analyze.py`). Coverage is binned by GC content, and each bin is scaled to remove the bias before calling:

```python
def gc_correct(cov_df):
    df = cov_df.copy()
    df["gc_bin"] = pd.qcut(df["gc_content"], q=20, duplicates="drop")
    gc_medians = df.groupby("gc_bin")["raw_coverage"].transform("median")
    global_median = df["raw_coverage"].median()
    df["corrected"] = df["raw_coverage"] * (global_median / gc_medians)
    return df["corrected"].values
```

## Files

```
generate_data.py   simulate per-bin coverage with planted CNVs and GC bias
analyze.py         GC correction, log2 ratio, rolling-median segmentation
plots.py           genome-wide log2-ratio track with gain/loss shading
```
