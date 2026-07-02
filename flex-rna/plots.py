#!/usr/bin/env python3
"""Plots for the Flex demo (writes assets/flex_qc.png)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    df = pd.read_csv(DATA / "cell_qc.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))

    # barcode-rank (knee) plot
    ranked = df.total_umi.sort_values(ascending=False).values
    ax[0].loglog(range(1, len(ranked) + 1), ranked, color=S.OUTLINE["blue"], linewidth=1.6)
    thr = df[df.called_cell].total_umi.min()
    ax[0].axhline(thr, ls="--", color=S.OUTLINE["pink"], lw=1.2, label="cell threshold")
    ax[0].set(title="Barcode rank (knee)", xlabel="barcode rank", ylabel="total UMI")
    ax[0].legend()

    # probes/cell per sample -- pastel filled boxes with fine outlines
    called = df[df.called_cell]
    samples = sorted(called.sample_bc.unique())
    data = [called[called.sample_bc == s].probes_detected.values for s in samples]
    keys = ["blue", "green", "peach", "lav"]
    bp = ax[1].boxplot(data, showfliers=False, patch_artist=True, widths=0.6,
                       medianprops=dict(color=S.INK, linewidth=1.2))
    for patch, k in zip(bp["boxes"], keys):
        patch.set(facecolor=S.PALETTE[k], edgecolor=S.OUTLINE[k], linewidth=1.1)
    for w in bp["whiskers"] + bp["caps"]:
        w.set(color=S.MUTED, linewidth=1.0)
    ax[1].set_xticks(range(1, len(samples) + 1)); ax[1].set_xticklabels(samples)
    ax[1].set(title="Probes detected / cell", ylabel="probes", xlabel="sample barcode")
    ax[1].grid(axis="x", visible=False)

    out = ASSETS / "flex_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
