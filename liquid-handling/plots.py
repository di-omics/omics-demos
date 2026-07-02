#!/usr/bin/env python3
"""Figures for the normalization demo: STARlet deck map + plate heatmaps.

Writes assets/normalization_qc.png (deck layout, per-well transfer volume, and
normalized output mass with flagged wells marked).
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.gridspec import GridSpec
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)
ROWS = list("ABCDEFGH")
INK, MUTED, TIP, PLATE, FLAG = "#14181D", "#5C6570", "#C77A3B", "#3B6EA5", "#C2367A"


def grid(df, value):
    return df.pivot(index="row", columns="col", values=value).reindex(ROWS).values


def deck_map(ax):
    ax.set_xlim(0, 22); ax.set_ylim(0, 6); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 22, 6, fc="#F1F3F5", ec="#D6DBE1"))
    ax.text(0.2, 5.6, "Hamilton STARlet deck", fontsize=9, color=MUTED, fontweight="bold")

    # tip carrier (rails 1): 3 racks
    ax.add_patch(FancyBboxPatch((0.6, 0.6), 4, 4.4, boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=TIP, lw=1.4))
    ax.text(2.6, 5.15, "tip carrier · rail 1", ha="center", fontsize=7.5, color=MUTED)
    for i, lab in enumerate(["tips: diluent", "tips: sample", "tips: pool"]):
        ax.add_patch(Rectangle((0.9, 0.9 + i * 1.35), 3.4, 1.1, fc=TIP, alpha=0.14, ec=TIP))
        ax.text(2.6, 1.45 + i * 1.35, lab, ha="center", va="center", fontsize=7.5, color=INK)

    # plate carrier (rails 10): 4 plates
    ax.add_patch(FancyBboxPatch((6, 0.6), 15, 4.4, boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=PLATE, lw=1.4))
    ax.text(13.5, 5.15, "plate carrier · rail 10", ha="center", fontsize=7.5, color=MUTED)
    plates = ["samples\n(source)", "diluent", "normalized\n(target)", "pool\n(8 fractions)"]
    for i, lab in enumerate(plates):
        ax.add_patch(Rectangle((6.4 + i * 3.6, 1.1), 3.1, 3.2, fc=PLATE, alpha=0.12, ec=PLATE))
        ax.text(6.4 + i * 3.6 + 1.55, 2.7, lab, ha="center", va="center", fontsize=7.5, color=INK)
        ax.text(6.4 + i * 3.6 + 1.55, 1.35, "96", ha="center", fontsize=6.5, color=MUTED)


def heat(ax, M, title, cmap, flags=None, unit=""):
    im = ax.imshow(M, cmap=cmap, aspect="auto")
    ax.set_xticks(range(12)); ax.set_xticklabels(range(1, 13), fontsize=7)
    ax.set_yticks(range(8)); ax.set_yticklabels(ROWS, fontsize=7)
    ax.set_title(title, fontsize=9)
    if flags is not None:
        for (r, c) in flags:
            ax.add_patch(Rectangle((c - 0.5, r - 0.5), 1, 1, fill=False, ec=FLAG, lw=1.8))
    cb = plt.colorbar(im, ax=ax, fraction=0.045, pad=0.03)
    cb.ax.tick_params(labelsize=7); cb.set_label(unit, fontsize=7)


def main():
    df = pd.read_csv(DATA / "normalization_report.tsv", sep="\t")
    flags = [(ROWS.index(r.row), r.col - 1) for _, r in df[df.status != "ok"].iterrows()]

    fig = plt.figure(figsize=(10, 6.2))
    gs = GridSpec(2, 2, height_ratios=[1.05, 1.25], hspace=0.42, wspace=0.28)
    deck_map(fig.add_subplot(gs[0, :]))
    heat(fig.add_subplot(gs[1, 0]), grid(df, "sample_ul"), "Sample volume per well",
         "Blues", unit="µL")
    heat(fig.add_subplot(gs[1, 1]), grid(df, "out_mass_ng"), "Normalized output mass (flagged outlined)",
         "viridis", flags=flags, unit="ng")
    fig.suptitle("liquid-handling - 96-well normalization on a Hamilton STARlet (PyLabRobot sim)",
                 fontsize=11, y=0.98)
    out = ASSETS / "normalization_qc.png"; fig.savefig(out, dpi=130, bbox_inches="tight",
                                                       facecolor="white")
    print("wrote", out)


if __name__ == "__main__":
    main()
