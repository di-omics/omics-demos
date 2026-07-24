#!/usr/bin/env python3
"""Plots for the scRNA-seq UMI deduplication demo.

Two panels: per-gene bar chart (three methods) + rarefaction curve showing
naive / directional / adjacency diverge.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t", index_col=0).head(10)
    sat = pd.read_csv(DATA / "saturation.tsv", sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))

    # per-gene bar chart: three methods
    x = range(len(counts))
    w = 0.25
    ax[0].bar([i - w for i in x], counts.naive_umi, width=w, label="naive",
              color=S.PALETTE["lav"], edgecolor=S.OUTLINE["lav"], linewidth=0.8)
    ax[0].bar(list(x), counts.directional_umi, width=w, label="directional",
              color=S.PALETTE["blue"], edgecolor=S.OUTLINE["blue"], linewidth=0.8)
    ax[0].bar([i + w for i in x], counts.adjacency_umi, width=w, label="adjacency",
              color=S.PALETTE["green"], edgecolor=S.OUTLINE["green"], linewidth=0.8)
    ax[0].set(title="Dedup methods compared (top genes)", ylabel="molecules",
              xticks=list(x))
    ax[0].set_xticklabels(counts.index, rotation=90, fontsize=7)
    ax[0].grid(axis="x", visible=False); ax[0].legend()

    # rarefaction / complexity curve
    ax[1].plot(sat.reads, sat.naive, marker="o", markersize=4,
               color=S.OUTLINE["lav"], markerfacecolor=S.PALETTE["lav"],
               markeredgecolor=S.OUTLINE["lav"], label="naive")
    ax[1].plot(sat.reads, sat.directional, marker="s", markersize=4,
               color=S.OUTLINE["blue"], markerfacecolor=S.PALETTE["blue"],
               markeredgecolor=S.OUTLINE["blue"], label="directional")
    ax[1].plot(sat.reads, sat.adjacency, marker="^", markersize=4,
               color=S.OUTLINE["green"], markerfacecolor=S.PALETTE["green"],
               markeredgecolor=S.OUTLINE["green"], label="adjacency")
    ax[1].set(title="Library complexity rarefaction",
              xlabel="reads sampled", ylabel="unique molecules")
    ax[1].legend()

    out = ASSETS / "scrnaseq_umi_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
