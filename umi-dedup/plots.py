#!/usr/bin/env python3
"""Plots for the UMI dedup demo (writes assets/umi_qc.png)."""
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t", index_col=0).head(10)
    sat = pd.read_csv(DATA / "saturation.tsv", sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    x = range(len(counts))
    ax[0].bar([i - 0.2 for i in x], counts.raw_reads, width=0.4, label="raw reads", color="#8a8a8a")
    ax[0].bar([i + 0.2 for i in x], counts.collapsed_umi, width=0.4, label="collapsed UMIs", color="#3b6ea5")
    ax[0].set(title="Raw reads vs deduped molecules", ylabel="count", xticks=list(x))
    ax[0].set_xticklabels(counts.index, rotation=90, fontsize=6)
    ax[0].legend(frameon=False)

    ax[1].plot(sat.reads, sat.molecules, marker="o", color="#a53b6e")
    ax[1].set(title="Library saturation", xlabel="reads sampled", ylabel="unique molecules")
    fig.tight_layout(); out = ASSETS / "umi_qc.png"; fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
