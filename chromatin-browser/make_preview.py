#!/usr/bin/env python3
"""Render a static preview of the chromatin browser for the README.

Mirrors the in-browser synthetic tracks (H3K27me3 broad domains, Pol II S5p
promoter peaks) and their called peaks. Open index.html for the interactive
version. All data synthetic.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from style.plot_style import set_style, finalize, PALETTE, OUTLINE, OUTLINE_WIDTH
set_style()

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from pathlib import Path

ASSETS = Path(__file__).parent / "assets"; ASSETS.mkdir(exist_ok=True)
LEN, N = 50000, 500
K27, POL2, PEAK = PALETTE["lavender"], PALETTE["mint"], PALETTE["rose"]
rng = np.random.default_rng(20260702)
bp = np.linspace(0, LEN, N)

k = np.zeros(N)
for amp, per, ph in [(1.0, 26000, 0.4), (0.6, 14000, 2.1), (0.35, 8000, 5.0)]:
    k += amp * (0.5 + 0.5 * np.sin(2 * np.pi * bp / per + ph))
k = np.clip(k / 1.95 + rng.normal(0, 0.02, N) - 0.28, 0, None); k /= k.max()

proms = [6200, 15800, 24500, 33900, 41200, 46800]
wid =   [900, 700, 1200, 650, 800, 600]
hei =   [0.95, 0.7, 1.0, 0.55, 0.8, 0.5]
p = 0.03 * rng.random(N)
for pr, w, h in zip(proms, wid, hei):
    p += h * np.exp(-((bp - pr) ** 2) / (2 * w * w))
p /= p.max()


def peaks(sig, thr):
    out, s = [], None
    for i, v in enumerate(sig):
        if v >= thr and s is None: s = i
        if (v < thr or i == N - 1) and s is not None:
            e = i - 1 if v < thr else i
            if e - s >= 1: out.append((bp[s], bp[e])); s = None
    return out


fig, axes = plt.subplots(2, 1, figsize=(9, 3.6), sharex=True)
for ax, sig, col, name, thr in [(axes[0], k, K27, "H3K27me3", 0.45),
                                (axes[1], p, POL2, "Pol II S5p", 0.30)]:
    ax.fill_between(bp / 1000, sig, color=col, alpha=0.25)
    ax.plot(bp / 1000, sig, color=col, lw=1.3)
    for a, b in peaks(sig, thr):
        ax.add_patch(FancyBboxPatch((a / 1000, 1.06), (b - a) / 1000, 0.12,
                     boxstyle="round,pad=0,rounding_size=0.3", fc=PEAK, ec="none",
                     alpha=0.85))
    ax.set_ylim(0, 1.28); ax.set_yticks([])
    ax.text(-0.5, 0.5, name, ha="right", va="center", fontsize=10, fontweight="bold")
    for s in ["top", "right", "left"]: ax.spines[s].set_visible(False)
axes[1].set_xlabel("position (kb)"); axes[1].set_xlim(0, 50)
fig.suptitle("chromatin-browser \u2013 synthetic CUT&Tag signal + called peaks",
             fontsize=11, x=0.55)
finalize(fig)
out = ASSETS / "preview.png"; fig.savefig(out, dpi=130)
print("wrote", out)
