#!/usr/bin/env python3
"""QC plots for the methylation-sequencing demo.

Four panels: per-CpG methylation histogram, spike-in controls,
CpG/CHG/CHH context breakdown, and M-bias (per-read-position).
"""
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
    fig, ax = plt.subplots(2, 2, figsize=(10, 7))

    # (0,0) Per-CpG methylation histogram
    for cond, f, key in [("10 ng", "cpg_10ng.tsv", "blue"),
                         ("0.1 ng", "cpg_0p1ng.tsv", "pink")]:
        d = pd.read_csv(DATA / f, sep="\t")
        ax[0, 0].hist(frac(d), bins=30, alpha=0.55, label=cond, density=True,
                       color=S.PALETTE[key], edgecolor=S.OUTLINE[key], linewidth=0.8)
    ax[0, 0].set(title="Per-CpG methylation", xlabel="methylation fraction",
                 ylabel="density")
    ax[0, 0].legend()

    # (0,1) Spike-in controls
    unmethylated = pd.read_csv(DATA / "control_unmethylated.tsv", sep="\t")
    methylated = pd.read_csv(DATA / "control_methylated.tsv", sep="\t")
    ce = 100 * (
        1
        - unmethylated.meth_reads.sum()
        / (unmethylated.meth_reads + unmethylated.unmeth_reads).sum()
    )
    pr = (
        100
        * methylated.meth_reads.sum()
        / (methylated.meth_reads + methylated.unmeth_reads).sum()
    )
    ax[0, 1].bar(["conversion", "protection"], [ce, pr],
                 color=[S.PALETTE["blue"], S.PALETTE["pink"]],
                 edgecolor=[S.OUTLINE["blue"], S.OUTLINE["pink"]], width=0.6)
    ax[0, 1].set(title="Spike-in controls", ylabel="%", ylim=(0, 114))
    ax[0, 1].grid(axis="x", visible=False)
    for i, v in enumerate([ce, pr]):
        ax[0, 1].text(i, v + 2, f"{v:.2f}%", ha="center", fontsize=8.5,
                      color=S.MUTED)

    # (1,0) CpG vs CHG vs CHH context breakdown
    contexts = ["CpG", "CHG", "CHH"]
    colors = [S.PALETTE["blue"], S.PALETTE["green"], S.PALETTE["lav"]]
    outlines = [S.OUTLINE["blue"], S.OUTLINE["green"], S.OUTLINE["lav"]]
    x_pos = range(len(contexts))
    bar_w = 0.35
    for i, (cond, f) in enumerate([("10 ng", "meth_10ng.tsv"),
                                    ("0.1 ng", "meth_0p1ng.tsv")]):
        d = pd.read_csv(DATA / f, sep="\t")
        rates = []
        for ctx in contexts:
            sub = d[d.context == ctx]
            m, u = sub.meth_reads.sum(), sub.unmeth_reads.sum()
            rates.append(100 * m / (m + u))
        offset = -bar_w / 2 + i * bar_w
        bars = ax[1, 0].bar([p + offset for p in x_pos], rates, width=bar_w,
                            label=cond, color=colors, alpha=0.7 if i == 0 else 0.5,
                            edgecolor=outlines, linewidth=0.8)
    ax[1, 0].set(title="Methylation by context", ylabel="% methylated",
                 xticks=list(x_pos))
    ax[1, 0].set_xticklabels(contexts)
    ax[1, 0].legend()
    ax[1, 0].grid(axis="x", visible=False)

    # (1,1) M-bias plot
    for cond, f, key in [("10 ng", "mbias_10ng.tsv", "blue"),
                         ("0.1 ng", "mbias_0p1ng.tsv", "pink")]:
        mb = pd.read_csv(DATA / f, sep="\t")
        ax[1, 1].plot(mb.read_position, mb.meth_rate * 100,
                      color=S.OUTLINE[key], linewidth=1.2, label=cond)
        ax[1, 1].fill_between(mb.read_position, mb.meth_rate * 100,
                              color=S.PALETTE[key], alpha=0.25)
    ax[1, 1].set(title="M-bias (per-read position)", xlabel="read position (bp)",
                 ylabel="% methylated")
    ax[1, 1].legend()

    out = ASSETS / "qc_plots.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
