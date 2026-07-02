#!/usr/bin/env python3
"""Plots for the UMI dedup demo (writes assets/umi_qc.png)."""
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

    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
    x = range(len(counts))
    ax[0].bar([i - 0.2 for i in x], counts.raw_reads, width=0.4, label="raw reads",
              color=S.PALETTE["lav"], edgecolor=S.OUTLINE["lav"], linewidth=0.9)
    ax[0].bar([i + 0.2 for i in x], counts.collapsed_umi, width=0.4, label="collapsed UMIs",
              color=S.PALETTE["blue"], edgecolor=S.OUTLINE["blue"], linewidth=0.9)
    ax[0].set(title="Raw reads vs deduped molecules", ylabel="count", xticks=list(x))
    ax[0].set_xticklabels(counts.index, rotation=90, fontsize=7)
    ax[0].grid(axis="x", visible=False); ax[0].legend()

    ax[1].plot(sat.reads, sat.molecules, marker="o", color=S.OUTLINE["pink"],
               markerfacecolor=S.PALETTE["pink"], markeredgecolor=S.OUTLINE["pink"], markersize=6)
    ax[1].set(title="Library saturation", xlabel="reads sampled", ylabel="unique molecules")

    out = ASSETS / "umi_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
