#!/usr/bin/env python3
"""QC plots for the EM-seq demo (writes assets/qc_plots.png)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def frac(d):
    return d.meth_reads / (d.meth_reads + d.unmeth_reads)


def main():
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.7))

    for cond, f, key in [("10 ng", "cpg_10ng.tsv", "blue"), ("0.1 ng", "cpg_0p1ng.tsv", "pink")]:
        d = pd.read_csv(DATA / f, sep="\t")
        ax[0].hist(frac(d), bins=30, alpha=0.55, label=cond, density=True,
                   color=S.PALETTE[key], edgecolor=S.OUTLINE[key], linewidth=0.8)
    ax[0].set(title="Per-CpG methylation", xlabel="methylation fraction", ylabel="density")
    ax[0].legend()

    lam = pd.read_csv(DATA / "control_lambda.tsv", sep="\t")
    puc = pd.read_csv(DATA / "control_puc19.tsv", sep="\t")
    ce = 100 * (1 - lam.meth_reads.sum() / (lam.meth_reads + lam.unmeth_reads).sum())
    pr = 100 * puc.meth_reads.sum() / (puc.meth_reads + puc.unmeth_reads).sum()
    ax[1].bar(["conversion\n(lambda)", "protection\n(pUC19)"], [ce, pr],
              color=[S.PALETTE["blue"], S.PALETTE["pink"]],
              edgecolor=[S.OUTLINE["blue"], S.OUTLINE["pink"]], width=0.6)
    ax[1].set(title="Spike-in controls", ylabel="%", ylim=(0, 114))
    ax[1].grid(axis="x", visible=False)
    for i, v in enumerate([ce, pr]):
        ax[1].text(i, v + 2, f"{v:.2f}%", ha="center", fontsize=8.5, color=S.MUTED)

    out = ASSETS / "qc_plots.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
