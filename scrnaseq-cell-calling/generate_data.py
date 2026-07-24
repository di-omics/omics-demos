#!/usr/bin/env python3
"""Synthesize a probe-based multiplexed scRNA-seq count matrix.

The generic design uses barcoded cells, probe-level UMI counts, and samples
multiplexed by sample barcode. It adds ambient background and empty barcodes for
downstream cell-calling and per-sample QC.

All data is synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(11)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

N_PROBES = 200
SAMPLES = ["BC01", "BC02", "BC03", "BC04"]     # multiplexed sample barcodes
CELLS_PER_SAMPLE = 300
EMPTY_FRAC = 0.25                               # background/empty barcodes


if __name__ == "__main__":
    probes = [f"probe_{i:03d}" for i in range(N_PROBES)]
    # each sample expresses a slightly different probe program
    programs = {s: RNG.gamma(1.5, 1.0, N_PROBES) for s in SAMPLES}

    cell_rows, meta = [], []
    cid = 0
    for s in SAMPLES:
        for _ in range(CELLS_PER_SAMPLE):
            depth = int(RNG.lognormal(7.6, 0.5))            # UMIs/cell
            lam = programs[s] / programs[s].sum() * depth
            counts = RNG.poisson(lam)
            cell_rows.append(counts); meta.append((f"cell_{cid}", s, "cell")); cid += 1
    # empty/background barcodes: low depth, ambient (avg program)
    n_empty = int(len(cell_rows) * EMPTY_FRAC)
    ambient = np.mean([programs[s] for s in SAMPLES], axis=0)
    for _ in range(n_empty):
        depth = int(RNG.lognormal(4.5, 0.6))
        lam = ambient / ambient.sum() * depth
        cell_rows.append(RNG.poisson(lam))
        meta.append((f"cell_{cid}", RNG.choice(SAMPLES), "empty")); cid += 1

    mat = pd.DataFrame(cell_rows, columns=probes)
    mat.insert(0, "barcode", [m[0] for m in meta])
    mat.to_csv(OUT / "counts.tsv", sep="\t", index=False)
    pd.DataFrame(meta, columns=["barcode", "sample_bc", "truth"]).to_csv(
        OUT / "barcodes.tsv", sep="\t", index=False)
    print(f"wrote {mat.shape[0]} barcodes x {N_PROBES} probes -> {OUT/'counts.tsv'}")
