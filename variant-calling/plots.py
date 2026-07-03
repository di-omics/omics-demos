#!/usr/bin/env python3
"""Plots for the variant-calling demo (writes assets/variant_qc.png).

Two panels: sensitivity curve by AF bin + true-vs-observed allele-fraction scatter.
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
    bins = pd.read_csv(DATA / "sensitivity_by_bin.tsv", sep="\t")
    merged = pd.read_csv(DATA / "truth_vs_calls.tsv", sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(11, 4))

    # --- left: sensitivity curve by AF bin ---
    x = range(len(bins))
    ax[0].plot(list(x), bins.sensitivity, marker="o", markersize=6,
               color=S.OUTLINE["pink"], markerfacecolor=S.PALETTE["pink"],
               markeredgecolor=S.OUTLINE["pink"], linewidth=1.6)
    ax[0].set(title="Sensitivity by allele-fraction bin",
              xlabel="allele fraction bin", ylabel="recall / sensitivity",
              ylim=(-0.05, 1.05), xticks=list(x))
    short = bins.af_bin.str.replace("0.", ".", regex=False)
    ax[0].set_xticklabels(short, rotation=30, ha="right", fontsize=7.5)
    ax[0].grid(axis="x", visible=False)

    # --- right: true AF vs observed AF scatter ---
    tp = merged[merged.called].copy()
    ax[1].plot([0, 1], [0, 1], linestyle="--", linewidth=1, color=S.MUTED,
               label="y = x", zorder=1)
    ax[1].scatter(tp.true_af, tp.obs_af, s=36, alpha=0.8, zorder=2,
                  color=S.PALETTE["green"], edgecolors=S.OUTLINE["green"],
                  linewidths=0.8, label="true positive")
    ax[1].set(title="True vs observed allele fraction",
              xlabel="true AF", ylabel="observed AF",
              xlim=(-0.02, 1.05), ylim=(-0.02, 1.05))
    ax[1].legend(fontsize=8)

    out = ASSETS / "variant_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
