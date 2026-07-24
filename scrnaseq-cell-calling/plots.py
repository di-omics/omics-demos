#!/usr/bin/env python3
"""Plots for the scRNA-seq cell-calling demo.

Three panels: barcode-rank (knee) with curvature-detected threshold,
probes/cell per sample, and ambient RNA profile.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    df = pd.read_csv(DATA / "cell_qc.tsv", sep="\t")
    ambient = pd.read_csv(DATA / "ambient_profile.tsv", sep="\t")

    fig, ax = plt.subplots(1, 3, figsize=(13, 4))

    # barcode-rank (knee) plot
    ranked = df.total_umi.sort_values(ascending=False).values
    ax[0].loglog(range(1, len(ranked) + 1), ranked, color=S.OUTLINE["blue"],
                 linewidth=1.6)
    thr = df[df.called_cell].total_umi.min()
    ax[0].axhline(thr, ls="--", color=S.OUTLINE["pink"], lw=1.2,
                  label=f"knee threshold ({thr:.0f})")
    n_cells = df.called_cell.sum()
    ax[0].axvline(n_cells, ls=":", color=S.MUTED, lw=0.8, alpha=0.6)
    ax[0].set(title="Barcode rank (max-curvature knee)",
              xlabel="barcode rank", ylabel="total UMI")
    ax[0].legend(fontsize=7)

    # UMI/cell per sample
    called = df[df.called_cell]
    samples = sorted(called.sample_bc.unique())
    data = [called[called.sample_bc == s].total_umi.values for s in samples]
    keys = ["blue", "green", "peach", "lav"]
    bp = ax[1].boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
                       medianprops=dict(color=S.INK, linewidth=1.2))
    for patch, k in zip(bp["boxes"], keys):
        patch.set(facecolor=S.PALETTE[k], edgecolor=S.OUTLINE[k], linewidth=1.1)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=S.MUTED, linewidth=1.0)
    ax[1].set_xticks(range(1, len(samples) + 1)); ax[1].set_xticklabels(samples)
    ax[1].set(title="UMI / cell", ylabel="total UMI",
              xlabel="sample barcode")
    ax[1].grid(axis="x", visible=False)

    # ambient RNA profile (top 20 probes)
    top = ambient.nlargest(20, "ambient_fraction")
    ax[2].barh(range(len(top)), top.ambient_fraction.values,
               color=S.PALETTE["teal"], edgecolor=S.OUTLINE["teal"], linewidth=0.8)
    ax[2].set_yticks(range(len(top)))
    ax[2].set_yticklabels(top.probe.values, fontsize=6.5)
    ax[2].set(title="Ambient RNA profile (top 20 probes)",
              xlabel="fraction of ambient UMI")
    ax[2].invert_yaxis()
    ax[2].grid(axis="y", visible=False)

    out = ASSETS / "scrnaseq_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
