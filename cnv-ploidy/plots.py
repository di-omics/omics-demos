#!/usr/bin/env python3
"""Plots for the CNV ploidy demo (writes assets/cnv_qc.png).

Single wide panel: genome-wide log2-ratio track with raw dots, smoothed line,
chromosome boundaries, and shaded gain/loss regions.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    cov = pd.read_csv(DATA / "normalized.tsv", sep="\t")
    called = pd.read_csv(DATA / "called_segments.tsv", sep="\t")

    # build genome-wide x coordinate (concatenated chromosomes)
    chroms = cov["chrom"].unique()
    offsets = {}
    x_global = np.zeros(len(cov), dtype=int)
    cum = 0
    boundaries = []
    chrom_mids = {}
    for chrom in chroms:
        mask = cov["chrom"] == chrom
        n_bins = mask.sum()
        offsets[chrom] = cum
        x_global[mask.values] = cov.loc[mask, "bin"].values + cum
        chrom_mids[chrom] = cum + n_bins // 2
        cum += n_bins
        boundaries.append(cum)

    cov["x"] = x_global

    fig, ax = plt.subplots(figsize=(12, 4))

    # raw log2 ratio as small dots
    ax.scatter(cov["x"], cov["log2"], s=3, alpha=0.5,
               color=S.PALETTE["blue"], linewidths=0, zorder=2)

    # smoothed line
    ax.plot(cov["x"], cov["smoothed"], linewidth=1.8,
            color=S.OUTLINE["blue"], zorder=3)

    # shade called regions
    for _, seg in called.iterrows():
        x_start = seg.start_bin + offsets[seg.chrom]
        x_end = seg.end_bin + offsets[seg.chrom]
        color = S.PALETTE["pink"] if seg["type"] == "gain" else S.PALETTE["lav"]
        ax.axvspan(x_start, x_end, alpha=0.35, color=color, zorder=1,
                   label=seg["type"])

    # chromosome boundaries
    for b in boundaries[:-1]:
        ax.axvline(b, color=S.HAIR, linewidth=0.8, zorder=1)

    # chromosome labels
    ax.set_xticks([chrom_mids[c] for c in chroms])
    ax.set_xticklabels(chroms)

    # reference lines
    ax.axhline(0, color=S.HAIR, linewidth=0.8, linestyle="--")

    ax.set(title="Genome-wide CNV profile", xlabel="chromosome",
           ylabel="log2 ratio")
    ax.set_xlim(0, cum)

    # deduplicate legend entries
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    ax.legend(by_label.values(), by_label.keys(), loc="upper right")

    out = ASSETS / "cnv_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
