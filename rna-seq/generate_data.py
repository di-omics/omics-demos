#!/usr/bin/env python3
"""Synthesize a bulk RNA-seq count matrix: two conditions, a planted DE gene set.

Gene baselines are drawn log-normally; a subset gets a real log2 fold-change in the
treated arm. Counts are negative-binomial (gamma-Poisson) with per-sample library
size variation, so the demo has to normalize before testing. All data is synthetic.
"""
import numpy as np, pandas as pd
from pathlib import Path

RNG = np.random.default_rng(13)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

N_GENES, N_DE, N_REP, DISP = 2000, 120, 5, 0.12

if __name__ == "__main__":
    genes = [f"gene_{i:04d}" for i in range(N_GENES)]
    base = RNG.lognormal(3.6, 1.2, N_GENES)               # per-gene baseline mean

    log2fc = np.zeros(N_GENES)
    de_idx = RNG.choice(N_GENES, N_DE, replace=False)
    fc = RNG.normal(0, 2.4, N_DE)
    fc[np.abs(fc) < 1.3] += np.sign(fc[np.abs(fc) < 1.3] + 1e-9) * 1.3   # push clear of |1|
    log2fc[de_idx] = fc

    samples, cond = [], []
    for c in ["ctrl", "treat"]:
        for r in range(N_REP):
            samples.append(f"{c}_{r+1}"); cond.append(c)
    lib = RNG.uniform(0.8, 1.25, len(samples))            # library-size factors

    mat = np.zeros((N_GENES, len(samples)), dtype=int)
    for s, (c, lf) in enumerate(zip(cond, lib)):
        eff = base * (2.0 ** (log2fc if c == "treat" else 0)) * lf
        rate = RNG.gamma(shape=1.0 / DISP, scale=eff * DISP)   # gamma-Poisson = NB
        mat[:, s] = RNG.poisson(rate)

    pd.DataFrame(mat, index=genes, columns=samples).rename_axis("gene").to_csv(
        OUT / "counts.tsv", sep="\t")
    pd.DataFrame({"sample": samples, "condition": cond}).to_csv(
        OUT / "samples.tsv", sep="\t", index=False)
    pd.DataFrame({"gene": [genes[i] for i in de_idx], "log2fc_true": log2fc[de_idx].round(3)}
                 ).to_csv(OUT / "truth_de.tsv", sep="\t", index=False)
    print(f"wrote {N_GENES} genes x {len(samples)} samples, {N_DE} planted DE -> {OUT/'counts.tsv'}")
