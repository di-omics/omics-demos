#!/usr/bin/env python3
"""Synthesize per-well DNA concentrations for a full 96-well plate.

The spread is deliberately realistic: most wells land in a normal range, a few are
too dilute to reach the target mass at a fixed final volume, and a few are so
concentrated they'd need a sub-microliter transfer. The normalization step has to
handle all three. All data is synthetic.
"""
import numpy as np, pandas as pd
from pathlib import Path

RNG = np.random.default_rng(5)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

ROWS = list("ABCDEFGH")
COLS = list(range(1, 13))


def well_order():
    # column-major: A1..H1, A2..H2, ... (matches an 8-channel head moving by column)
    return [f"{r}{c}" for c in COLS for r in ROWS]


if __name__ == "__main__":
    wells = well_order()
    n = len(wells)                                  # 96
    conc = RNG.uniform(6, 45, n)                    # bulk of the plate: normal range
    dilute_idx = RNG.choice(n, 5, replace=False)    # a few too-dilute wells
    conc[dilute_idx] = RNG.uniform(1.4, 2.2, 5)
    pool = [i for i in range(n) if i not in dilute_idx]
    conc_idx = RNG.choice(pool, 6, replace=False)
    conc[conc_idx] = RNG.uniform(90, 115, 6)        # a few very concentrated wells
    conc = conc.round(1)

    df = pd.DataFrame({
        "well": wells,
        "row": [w[0] for w in wells],
        "col": [int(w[1:]) for w in wells],
        "conc_ng_ul": conc,
    })
    df.to_csv(OUT / "concentrations.tsv", sep="\t", index=False)
    print(f"wrote {n}-well plate -> {OUT/'concentrations.tsv'}")
    print(f"range {conc.min()}-{conc.max()} ng/uL "
          f"({int((conc < 2.1).sum())} dilute, {int((conc > 89).sum())} concentrated)")
