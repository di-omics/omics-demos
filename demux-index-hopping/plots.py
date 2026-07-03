#!/usr/bin/env python3
"""Plots for the demux index-hopping demo (writes assets/demux_qc.png).

Two panels: i7 x i5 heatmap + per-sample assigned-read bar chart.
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
    mat = pd.read_csv(DATA / "hopping_matrix.csv", index_col=0)
    per_sample = pd.read_csv(DATA / "per_sample.csv")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))

    # --- left: i7 x i5 heatmap ---
    vals = mat.values.astype(float)
    im = ax[0].imshow(vals, cmap=S.cmap("blue"), aspect="equal")
    n = len(mat)
    for i in range(n):
        for j in range(n):
            v = int(vals[i, j])
            if v > 0:
                ax[0].text(j, i, str(v), ha="center", va="center", fontsize=7,
                           color="white" if i == j else S.OUTLINE["blue"])
    ax[0].set_xticks(range(n)); ax[0].set_yticks(range(n))
    ax[0].set_xticklabels(mat.columns, fontsize=8)
    ax[0].set_yticklabels(mat.index, fontsize=8)
    ax[0].set_xlabel("i5 sample"); ax[0].set_ylabel("i7 sample")
    ax[0].set_title("Index-pair count matrix (i7 x i5)")
    fig.colorbar(im, ax=ax[0], fraction=0.046, pad=0.04)

    # --- right: per-sample bar chart ---
    x = range(len(per_sample))
    ax[1].bar(x, per_sample.assigned_reads,
              color=S.PALETTE["blue"], edgecolor=S.OUTLINE["blue"], linewidth=0.8)
    ax[1].set_xticks(list(x))
    ax[1].set_xticklabels(per_sample["sample"], fontsize=8)
    ax[1].set(title="Assigned reads per sample", ylabel="reads", xlabel="sample")
    ax[1].grid(axis="x", visible=False)

    out = ASSETS / "demux_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
