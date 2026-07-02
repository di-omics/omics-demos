#!/usr/bin/env python3
"""Plots for the Flex demo (writes assets/flex_qc.png)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from style.plot_style import set_style, finalize, PALETTE, OUTLINE, OUTLINE_WIDTH
set_style()

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    df = pd.read_csv(DATA / "cell_qc.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))

    # barcode-rank (knee) plot
    ranked = df.total_umi.sort_values(ascending=False).values
    ax[0].loglog(range(1, len(ranked) + 1), ranked, color=PALETTE["blue"])
    thr = df[df.called_cell].total_umi.min()
    ax[0].axhline(thr, ls="--", color=PALETTE["rose"], lw=1, label="cell threshold")
    ax[0].set(title="Barcode rank (knee)", xlabel="barcode rank", ylabel="total UMI")
    ax[0].legend()

    # probes/cell per sample
    called = df[df.called_cell]
    samples = sorted(called.sample_bc.unique())
    data = [called[called.sample_bc == s].probes_detected.values for s in samples]
    bp = ax[1].boxplot(data, showfliers=False, patch_artist=True,
                       medianprops=dict(color=PALETTE["ink"], linewidth=1))
    colors = [PALETTE["blue"], PALETTE["rose"], PALETTE["mint"], PALETTE["lavender"]]
    for patch, c in zip(bp["boxes"], colors):
        patch.set_facecolor(c)
        patch.set_edgecolor(OUTLINE)
        patch.set_linewidth(OUTLINE_WIDTH)
    ax[1].set_xticks(range(1, len(samples) + 1)); ax[1].set_xticklabels(samples)
    ax[1].set(title="Probes detected / cell", ylabel="probes", xlabel="sample barcode")
    finalize(fig)
    out = ASSETS / "flex_qc.png"; fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
