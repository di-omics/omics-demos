#!/usr/bin/env python3
"""CNV calling from binned coverage data.

Steps:
  1. Normalize raw coverage by genome-wide median.
  2. Correct for GC bias (rolling median of coverage vs GC).
  3. Compute log2 ratio.
  4. Smooth with a rolling median (window = 11 bins, per chromosome).
  5. Call gains (smoothed log2 > 0.3) and losses (smoothed log2 < -0.3).
  6. Compare called segments to planted truth.
"""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"

SMOOTH_WINDOW = 11
GAIN_THRESH = 0.3
LOSS_THRESH = -0.3


def gc_correct(cov_df):
    """Simple GC-bias correction: bin GC into quantiles, compute median
    coverage per bin, divide out."""
    df = cov_df.copy()
    df["gc_bin"] = pd.qcut(df["gc_content"], q=20, duplicates="drop")
    gc_medians = df.groupby("gc_bin")["raw_coverage"].transform("median")
    global_median = df["raw_coverage"].median()
    df["corrected"] = df["raw_coverage"] * (global_median / gc_medians)
    return df["corrected"].values


def call_segments(chrom_series, threshold_gain, threshold_loss):
    """Return list of (start_bin, end_bin, type) from a boolean mask."""
    segments = []
    in_event = False
    for idx, (b, val) in enumerate(chrom_series.items()):
        if val > threshold_gain:
            etype = "gain"
        elif val < threshold_loss:
            etype = "loss"
        else:
            etype = None

        if etype and not in_event:
            start = b
            cur_type = etype
            in_event = True
        elif etype and in_event and etype != cur_type:
            segments.append((start, b, cur_type))
            start = b
            cur_type = etype
        elif not etype and in_event:
            segments.append((start, b, cur_type))
            in_event = False
    if in_event:
        segments.append((start, b + 1, cur_type))
    return segments


def main():
    cov = pd.read_csv(DATA / "coverage.tsv", sep="\t")
    truth = pd.read_csv(DATA / "truth.tsv", sep="\t")

    # GC correction
    cov["corrected"] = gc_correct(cov)

    # normalize by median
    median_cov = cov["corrected"].median()
    cov["ratio"] = cov["corrected"] / median_cov

    # log2 ratio
    cov["log2"] = np.log2(cov["ratio"].clip(lower=1e-6))

    # smooth per chromosome
    cov["smoothed"] = (
        cov.groupby("chrom")["log2"]
        .transform(lambda s: s.rolling(SMOOTH_WINDOW, center=True, min_periods=1).median())
    )

    # call segments
    called = []
    for chrom, g in cov.groupby("chrom"):
        segs = call_segments(
            g.set_index("bin")["smoothed"], GAIN_THRESH, LOSS_THRESH
        )
        for s, e, t in segs:
            called.append((chrom, s, e, t))

    called_df = pd.DataFrame(called, columns=["chrom", "start_bin", "end_bin", "type"])

    # compare to truth
    print("\n=== CNV ploidy QC (synthetic) ===")
    print(f"Total bins:       {len(cov)}")
    print(f"Median coverage:  {median_cov:.1f}")
    print(f"Segments called:  {len(called_df)}")
    print()

    for _, ev in truth.iterrows():
        match = called_df[
            (called_df.chrom == ev.chrom) &
            (called_df.type == ev["type"]) &
            (called_df.start_bin <= ev.end_bin) &
            (called_df.end_bin >= ev.start_bin)
        ]
        status = "RECOVERED" if len(match) > 0 else "MISSED"
        print(f"  {ev.chrom}:{ev.start_bin}-{ev.end_bin} {ev['type']} "
              f"(ratio {ev.copy_ratio}) -> {status}")

    print()
    print("Called segments:")
    print(called_df.to_string(index=False))

    # save results
    cov.to_csv(DATA / "normalized.tsv", sep="\t", index=False)
    called_df.to_csv(DATA / "called_segments.tsv", sep="\t", index=False)
    print(f"\nwrote {DATA / 'normalized.tsv'} and {DATA / 'called_segments.tsv'}")


if __name__ == "__main__":
    main()
