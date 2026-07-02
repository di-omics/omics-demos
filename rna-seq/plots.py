#!/usr/bin/env python3
"""Plots for the RNA-seq demo: PCA scatter + volcano (writes assets/rna_qc.png)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    pca = pd.read_csv(DATA / "pca.tsv", sep="\t")
    de = pd.read_csv(DATA / "de_results.tsv", sep="\t")

    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.8))

    # PCA scatter
    colors = {"control": "#3b6ea5", "treatment": "#a53b6e"}
    for cond, grp in pca.groupby("condition"):
        ax[0].scatter(grp.PC1, grp.PC2, c=colors[cond], label=cond, s=60, edgecolors="w", lw=0.5)
    ax[0].set(title="PCA (log2-CPM)", xlabel="PC1", ylabel="PC2")
    ax[0].legend(frameon=False)

    # Volcano
    neg_log10p = -np.log10(de.pvalue.clip(lower=1e-300))
    sig = de.padj < 0.05
    ax[1].scatter(de.log2FC[~sig], neg_log10p[~sig], c="#999999", s=8, alpha=0.5, label="NS")
    ax[1].scatter(de.log2FC[sig], neg_log10p[sig], c="#d9534f", s=12, alpha=0.7, label="padj < 0.05")
    ax[1].set(title="Volcano plot", xlabel="log2 fold-change", ylabel="-log10(p)")
    ax[1].legend(frameon=False, markerscale=1.5)

    fig.tight_layout()
    out = ASSETS / "rna_qc.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
