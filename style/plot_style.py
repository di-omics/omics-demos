"""Shared plot style for omics-demos.

Registers bundled Manrope TTF fonts, sets a baby-pastel palette with fine
outlines, and provides a finalize() helper that prevents text overlap.

Usage in any demo's plots.py:
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
    from style.plot_style import set_style, finalize, PALETTE
    set_style()
    # ... build figure ...
    finalize(fig)
"""
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

# ---------------------------------------------------------------------------
# Font registration: bundle Manrope .ttf (converted from @fontsource woff2
# with name-table fixes so matplotlib sees "Manrope" instead of "ExtraLight")
# ---------------------------------------------------------------------------
_FONT_DIR = Path(__file__).parent / "fonts"

def _register_fonts():
    for ttf in sorted(_FONT_DIR.glob("*.ttf")):
        fm.fontManager.addfont(str(ttf))
    resolved = fm.findfont("Manrope", fallback_to_default=False)
    if "Manrope" not in resolved:
        raise RuntimeError(
            f"Manrope font not resolved (got {resolved}). "
            "Try: rm ~/.cache/matplotlib/*.json"
        )

_register_fonts()

# ---------------------------------------------------------------------------
# Baby-pastel palette with named entries
# ---------------------------------------------------------------------------
PALETTE = {
    "blue":     "#A8C8E8",
    "rose":     "#E8A8C0",
    "mint":     "#A8E8D0",
    "lavender": "#C8B8E8",
    "peach":    "#E8C8A8",
    "sky":      "#A8D8E8",
    "lilac":    "#D8A8E8",
    "sage":     "#C0D8B0",
    "coral":    "#E8B0A8",
    "cream":    "#E8E0C8",
    "grey":     "#C8C8C8",
    "ink":      "#14181D",
    "muted":    "#5C6570",
    "bg":       "#FAFBFC",
}

# Ordered list for cycling through series
PALETTE_CYCLE = [
    PALETTE["blue"], PALETTE["rose"], PALETTE["mint"], PALETTE["lavender"],
    PALETTE["peach"], PALETTE["sky"], PALETTE["lilac"], PALETTE["sage"],
]

# Fine outline color for patches/spines
OUTLINE = "#6B7280"
OUTLINE_WIDTH = 0.7

# ---------------------------------------------------------------------------
# set_style() - call once before creating figures
# ---------------------------------------------------------------------------
def set_style():
    plt.rcParams.update({
        "font.family": "Manrope",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.titleweight": 600,
        "axes.labelsize": 9.5,
        "axes.linewidth": OUTLINE_WIDTH,
        "axes.edgecolor": OUTLINE,
        "axes.facecolor": "white",
        "axes.prop_cycle": plt.cycler(color=PALETTE_CYCLE),
        "xtick.major.width": OUTLINE_WIDTH,
        "ytick.major.width": OUTLINE_WIDTH,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "figure.facecolor": PALETTE["bg"],
        "savefig.facecolor": PALETTE["bg"],
        "patch.linewidth": OUTLINE_WIDTH,
        "patch.edgecolor": OUTLINE,
        "lines.linewidth": 1.4,
    })

# ---------------------------------------------------------------------------
# finalize(fig) - run after building the figure; prevents text overlap
# ---------------------------------------------------------------------------
def finalize(fig):
    try:
        fig.set_layout_engine("constrained")
    except Exception:
        fig.tight_layout()
