#!/usr/bin/env python3
"""Synthesize a small targeted PCR sample sheet (one column, 8 samples).

Each sample gets an input amount and a unique i7/i5 index pair for the PCR2
indexing step. This drives a generic library-prep automation demo; the volumes
and layout in the protocol are illustrative, not a validated method.

All data is synthetic.
"""
import numpy as np, pandas as pd
from pathlib import Path

RNG = np.random.default_rng(5)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

ROWS = list("ABCDEFGH")
I7 = [f"i7_{n}" for n in ["N701", "N702", "N703", "N704", "N705", "N706", "N707", "N708"]]
I5 = [f"i5_{n}" for n in ["S507", "S508", "S510", "S511", "S513", "S515", "S516", "S517"]]

if __name__ == "__main__":
    df = pd.DataFrame({
        "well": [f"{r}1" for r in ROWS],
        "sample": [f"lib_{i+1:02d}" for i in range(8)],
        "input_ng": RNG.uniform(5, 25, 8).round(1),
        "i7_index": I7,
        "i5_index": I5,
    })
    df.to_csv(OUT / "sample_sheet.tsv", sep="\t", index=False)
    print(f"wrote {len(df)}-sample sheet (column 1) -> {OUT/'sample_sheet.tsv'}")
    print(df.to_string(index=False))
