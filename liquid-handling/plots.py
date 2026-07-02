#!/usr/bin/env python3
"""Figures for the liquid-handling demo: STARlet deck map + plate heatmaps.

Writes assets/normalization_qc.png.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from style.plot_style import set_style, finalize, PALETTE, OUTLINE, OUTLINE_WIDTH
set_style()

import numpy as np, pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.gridspec import GridSpec
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)
ROWS = list("ABCDEFGH")
FLAG = PALETTE["coral"]


def grid(df, value):
    return df.pivot(index="row", columns="col", values=value).reindex(ROWS).values


def deck_map(ax):
    ax.set_xlim(0, 22); ax.set_ylim(0, 6); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 22, 6, fc=PALETTE["bg"], ec=OUTLINE, lw=OUTLINE_WIDTH))
    ax.text(0.2, 5.6, "Hamilton STARlet deck", fontsize=9, color=PALETTE["muted"],
            fontweight="bold")

    # tip carrier (rails 1): 3 racks
    ax.add_patch(FancyBboxPatch((0.6, 0.6), 4, 4.4,
                 boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=PALETTE["peach"], lw=OUTLINE_WIDTH))
    ax.text(2.6, 5.15, "tip carrier", ha="center", fontsize=7.5, color=PALETTE["muted"])
    for i, lab in enumerate(["tips: diluent", "tips: sample", "tips: pool"]):
        ax.add_patch(Rectangle((0.9, 0.9 + i * 1.35), 3.4, 1.1,
                     fc=PALETTE["peach"], alpha=0.25, ec=OUTLINE, lw=OUTLINE_WIDTH))
        ax.text(2.6, 1.45 + i * 1.35, lab, ha="center", va="center",
                fontsize=7.5, color=PALETTE["ink"])

    # plate carrier (rails 10): 4 plates
    ax.add_patch(FancyBboxPatch((6, 0.6), 15, 4.4,
                 boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=PALETTE["blue"], lw=OUTLINE_WIDTH))
    ax.text(13.5, 5.15, "plate carrier", ha="center", fontsize=7.5, color=PALETTE["muted"])
    plates = ["samples\n(source)", "diluent", "normalized\n(target)", "pool\n(8 fractions)"]
    for i, lab in enumerate(plates):
        ax.add_patch(Rectangle((6.4 + i * 3.6, 1.1), 3.1, 3.2,
                     fc=PALETTE["blue"], alpha=0.18, ec=OUTLINE, lw=OUTLINE_WIDTH))
        ax.text(6.4 + i * 3.6 + 1.55, 2.7, lab, ha="center", va="center",
                fontsize=7.5, color=PALETTE["ink"])
        ax.text(6.4 + i * 3.6 + 1.55, 1.35, "96", ha="center", fontsize=6.5,
                color=PALETTE["muted"])


def heat(ax, M, title, cmap, flags=None, unit=""):
    im = ax.imshow(M, cmap=cmap, aspect="auto")
    ax.set_xticks(range(12)); ax.set_xticklabels(range(1, 13))
    ax.set_yticks(range(8)); ax.set_yticklabels(ROWS)
    ax.set_title(title)
    if flags is not None:
        for (r, c) in flags:
            ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False,
                         ec=FLAG, lw=1.6))
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=7); cb.set_label(unit, fontsize=7)


def main():
    df = pd.read_csv(DATA / "normalization_report.tsv", sep="\t")
    flags = [(ROWS.index(r.row), r.col - 1) for _, r in df[df.status != "ok"].iterrows()]

    fig = plt.figure(figsize=(10, 6.2))
    gs = GridSpec(2, 2, height_ratios=[1.05, 1.25], hspace=0.42, wspace=0.28)
    deck_map(fig.add_subplot(gs[0, :]))

    from matplotlib.colors import LinearSegmentedColormap
    blue_cmap = LinearSegmentedColormap.from_list("pblu", ["white", PALETTE["blue"]])
    mint_cmap = LinearSegmentedColormap.from_list("pmnt", [PALETTE["cream"], PALETTE["mint"]])

    heat(fig.add_subplot(gs[1, 0]), grid(df, "sample_ul"), "Sample volume per well",
         blue_cmap, unit="µL")
    heat(fig.add_subplot(gs[1, 1]), grid(df, "out_mass_ng"),
         "Normalized output mass (flagged outlined)",
         mint_cmap, flags=flags, unit="ng")
    fig.suptitle("liquid-handling \u2013 96-well normalization on a Hamilton STARlet",
                 fontsize=11, y=0.98)
    out = ASSETS / "normalization_qc.png"; fig.savefig(out, dpi=130, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
