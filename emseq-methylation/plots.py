#!/usr/bin/env python3
"""Minimal QC plots for the EM-seq demo (writes assets/qc_plots.png)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from style.plot_style import set_style, finalize, PALETTE, OUTLINE, OUTLINE_WIDTH
set_style()

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA = ROOT / "data"
ASSETS = ROOT / "assets"
ASSETS.mkdir(exist_ok=True)


def frac(d):
    return d.meth_reads / (d.meth_reads + d.unmeth_reads)


def main():
    fig, ax = plt.subplots(1, 2, figsize=(9, 3.4))

    for cond, f, c in [("10 ng", "cpg_10ng.tsv", PALETTE["blue"]),
                       ("0.1 ng", "cpg_0p1ng.tsv", PALETTE["rose"])]:
        d = pd.read_csv(DATA / f, sep="\t")
        ax[0].hist(frac(d), bins=30, alpha=0.7, label=cond, color=c,
                   density=True, edgecolor=OUTLINE, linewidth=OUTLINE_WIDTH)
    ax[0].set(title="Per-CpG methylation", xlabel="methylation fraction",
              ylabel="density")
    ax[0].legend()

    lam = pd.read_csv(DATA / "control_lambda.tsv", sep="\t")
    puc = pd.read_csv(DATA / "control_puc19.tsv", sep="\t")
    ce = 100 * (1 - lam.meth_reads.sum() / (lam.meth_reads + lam.unmeth_reads).sum())
    pr = 100 * puc.meth_reads.sum() / (puc.meth_reads + puc.unmeth_reads).sum()
    ax[1].bar(["conversion\n(lambda)", "protection\n(pUC19)"], [ce, pr],
              color=[PALETTE["blue"], PALETTE["rose"]],
              edgecolor=OUTLINE, linewidth=OUTLINE_WIDTH)
    ax[1].set(title="Spike-in controls", ylabel="%", ylim=(0, 112))
    for i, v in enumerate([ce, pr]):
        ax[1].text(i, v + 2, f"{v:.2f}%", ha="center", fontsize=9)

    finalize(fig)
    out = ASSETS / "qc_plots.png"
    fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
