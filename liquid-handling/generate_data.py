#!/usr/bin/env python3
"""Synthesize per-well DNA concentrations for the normalization demo (8 wells)."""
import numpy as np, pandas as pd
from pathlib import Path

RNG = np.random.default_rng(5)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)
WELLS = [f"{r}1" for r in "ABCDEFGH"]

if __name__ == "__main__":
    conc = RNG.uniform(4, 60, len(WELLS)).round(1)   # ng/uL, spread across a plate
    pd.DataFrame({"well": WELLS, "conc_ng_ul": conc}).to_csv(OUT / "concentrations.tsv", sep="\t", index=False)
    print("wrote", OUT / "concentrations.tsv")
    print(dict(zip(WELLS, conc)))
