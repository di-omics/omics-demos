#!/usr/bin/env python3
"""Plots for the RNA-seq demo: PCA of samples + a volcano of the DE test.
Writes assets/rnaseq_qc.png in the shared baby-pastel / Manrope style.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import numpy as np, pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)
COND = {"ctrl": "blue", "treat": "pink"}


def main():
    pca = pd.read_csv(DATA / "pca.tsv", sep="\t")
    res = pd.read_csv(DATA / "de_results.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.1))

    # PCA scatter
    for c, key in COND.items():
        sub = pca[pca.condition == c]
        ax[0].scatter(sub.PC1, sub.PC2, s=90, color=S.PALETTE[key],
                      edgecolor=S.OUTLINE[key], linewidth=1.2, label=c, zorder=3)
    ax[0].axhline(0, color=S.HAIR, lw=0.8); ax[0].axvline(0, color=S.HAIR, lw=0.8)
    ax[0].set(title="Sample PCA", xlabel="PC1", ylabel="PC2"); ax[0].legend()

    # volcano
    res = res.copy(); res["nlp"] = -np.log10(res.padj.clip(lower=1e-12))
    bg, hit = res[~res.called], res[res.called]
    ax[1].scatter(bg.log2fc, bg.nlp, s=10, color=S.GRID, edgecolor=S.HAIR, linewidth=0.3, zorder=2)
    ax[1].scatter(hit.log2fc, hit.nlp, s=18, color=S.PALETTE["pink"],
                  edgecolor=S.OUTLINE["pink"], linewidth=0.6, zorder=3, label="called DE")
    ax[1].axhline(-np.log10(0.05), ls="--", color=S.OUTLINE["teal"], lw=1.0)
    for v in (-1, 1):
        ax[1].axvline(v, ls="--", color=S.OUTLINE["teal"], lw=1.0)
    ax[1].set(title="Differential expression (volcano)", xlabel="log2 fold-change",
              ylabel="-log10 adjusted p"); ax[1].legend()

    fig.suptitle("rna-seq - bulk differential expression on synthetic counts",
                 fontsize=12, weight="semibold")
    out = ASSETS / "rnaseq_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
