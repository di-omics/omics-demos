#!/usr/bin/env python3
"""Synthesize per-bin coverage across multiple chromosomes with planted CNVs.

Models a shallow whole-genome sequencing experiment: ~100x baseline coverage
with Poisson noise, a few copy-number events (gains and losses), and mild
GC-bias (sinusoidal variation ~10%).

All data is synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(1)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

CHROMS = ["chr1", "chr2", "chr3", "chr4", "chr5"]
BINS_PER_CHROM = 200
BASELINE_COV = 100

# CNV events: (chrom, start_bin, end_bin, type, copy_ratio)
EVENTS = [
    ("chr2",  50,  80, "gain", 1.5),
    ("chr3", 100, 120, "loss", 0.5),
    ("chr5",  30,  50, "gain", 2.0),
]

if __name__ == "__main__":
    rows = []
    for chrom in CHROMS:
        for b in range(BINS_PER_CHROM):
            # GC content: smooth sinusoidal variation between 0.3 and 0.6
            gc = 0.45 + 0.15 * np.sin(2 * np.pi * b / BINS_PER_CHROM
                                       + RNG.uniform(0, 2 * np.pi))
            gc = np.clip(gc, 0.2, 0.8)

            # copy-ratio multiplier from planted events
            ratio = 1.0
            for ev_chrom, s, e, _, cr in EVENTS:
                if chrom == ev_chrom and s <= b < e:
                    ratio = cr
                    break

            # mild GC bias: ~10% sinusoidal modulation keyed to GC content
            gc_bias = 1.0 + 0.10 * np.sin(2 * np.pi * (gc - 0.3))

            expected = BASELINE_COV * ratio * gc_bias
            raw_cov = RNG.poisson(max(expected, 1))
            rows.append((chrom, b, raw_cov, round(gc, 4)))

    cov = pd.DataFrame(rows, columns=["chrom", "bin", "raw_coverage", "gc_content"])
    cov.to_csv(OUT / "coverage.tsv", sep="\t", index=False)
    print(f"wrote {len(cov)} bins -> {OUT / 'coverage.tsv'}")

    truth = pd.DataFrame(EVENTS,
                          columns=["chrom", "start_bin", "end_bin", "type", "copy_ratio"])
    truth.to_csv(OUT / "truth.tsv", sep="\t", index=False)
    print(f"wrote {len(truth)} planted CNV events -> {OUT / 'truth.tsv'}")
