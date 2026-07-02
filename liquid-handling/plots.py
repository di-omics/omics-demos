#!/usr/bin/env python3
"""Figures for the library-prep demo: Hamilton STAR deck map + per-step transfers.
Writes assets/libprep_qc.png in the shared baby-pastel / Manrope style.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import omics_style as S; S.apply()

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.gridspec import GridSpec

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)
PHASE_COLOR = {"PCR1 master mix": "blue", "SPRI cleanup": "teal", "PCR2 index MM": "pink"}


def box(ax, x, y, w, h, key, label, sub=None, fs=8.5):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.12",
                 fc=S.PALETTE[key], ec=S.OUTLINE[key], lw=1.1, alpha=0.9))
    ax.text(x + w / 2, y + h / 2 + (0.14 if sub else 0), label, ha="center", va="center",
            fontsize=fs, color=S.INK, weight="medium")
    if sub:
        ax.text(x + w / 2, y + h / 2 - 0.32, sub, ha="center", va="center",
                fontsize=6.8, color=S.MUTED)


def deck_map(ax):
    ax.set_xlim(0, 24); ax.set_ylim(0, 6); ax.axis("off")
    ax.add_patch(Rectangle((0, 0), 24, 6, fc="#F7F8FA", ec=S.HAIR, lw=1))
    ax.text(0.3, 5.55, "Hamilton STAR deck", fontsize=9.5, color=S.INK, weight="semibold")

    # tip carrier (2 racks)
    ax.add_patch(FancyBboxPatch((0.5, 0.5), 4.4, 4.4, boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=S.MUTED, lw=1))
    ax.text(2.7, 5.05, "tip carrier", ha="center", fontsize=7.2, color=S.MUTED)
    box(ax, 0.85, 2.75, 3.7, 1.75, "peach", "p50 tips")
    box(ax, 0.85, 0.85, 3.7, 1.75, "gold", "p300 tips")

    # plate carrier (5 plates)
    ax.add_patch(FancyBboxPatch((6, 0.5), 17.5, 4.4, boxstyle="round,pad=0.05,rounding_size=0.15",
                 fc="white", ec=S.MUTED, lw=1))
    ax.text(14.75, 5.05, "plate carrier", ha="center", fontsize=7.2, color=S.MUTED)
    plates = [("lav", "reagent", "PCR1 / PCR2 MM"), ("blue", "work", "samples + eluate"),
              ("teal", "magnet", "SPRI cleanup"), ("green", "reservoir", "beads / EtOH / EB"),
              ("pink", "waste", "spent liquid")]
    for i, (key, name, sub) in enumerate(plates):
        box(ax, 6.35 + i * 3.35, 1.05, 2.95, 3.3, key, name, sub=sub)


def transfers(ax):
    wl = pd.read_csv(DATA / "worklist.csv")
    wl = wl.iloc[::-1].reset_index(drop=True)   # first step on top
    y = range(len(wl))
    colors = [S.PALETTE[PHASE_COLOR[p]] for p in wl.phase]
    edges = [S.OUTLINE[PHASE_COLOR[p]] for p in wl.phase]
    ax.barh(list(y), wl.volume_ul, color=colors, edgecolor=edges, height=0.66)
    ax.set_yticks(list(y)); ax.set_yticklabels(wl.reagent, fontsize=8)
    ax.set_xlabel("volume per channel (µL)"); ax.set_xlim(0, wl.volume_ul.max() * 1.18)
    ax.set_title("Transfers by step (column-1, 8-channel)")
    ax.grid(axis="y", visible=False)
    for i, v in zip(y, wl.volume_ul):
        ax.text(v + wl.volume_ul.max() * 0.02, i, f"{v:g}", va="center", fontsize=7.5, color=S.MUTED)
    # phase legend
    handles = [plt.Line2D([0], [0], marker="s", ls="", markersize=8,
               mfc=S.PALETTE[c], mec=S.OUTLINE[c]) for c in PHASE_COLOR.values()]
    ax.legend(handles, list(PHASE_COLOR.keys()), loc="lower right", fontsize=8)


def main():
    fig = plt.figure(figsize=(10, 7))
    gs = GridSpec(2, 1, height_ratios=[1, 1.35], figure=fig)
    deck_map(fig.add_subplot(gs[0]))
    transfers(fig.add_subplot(gs[1]))
    fig.suptitle("liquid-handling - amplicon-seq library prep on a Hamilton STAR (PyLabRobot sim)",
                 fontsize=12, weight="semibold")
    out = ASSETS / "libprep_qc.png"; fig.savefig(out)
    print("wrote", out)


if __name__ == "__main__":
    main()
