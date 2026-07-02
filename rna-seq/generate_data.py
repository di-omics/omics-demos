#!/usr/bin/env python3
"""Synthesize a gene x sample count matrix for bulk/single-blastocyst RNA-seq.

Two conditions (control vs treatment), 4 replicates each. A planted set of DE
genes is spiked in so downstream analysis can measure recovery. All data is
synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(42)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

N_GENES = 2000
N_REPS = 4
CONDITIONS = ["control", "treatment"]
N_DE = 80  # planted differentially-expressed genes (40 up, 40 down)
DE_FOLD = 3.0  # fold-change for planted DE genes

if __name__ == "__main__":
    genes = [f"gene_{i:04d}" for i in range(N_GENES)]
    samples = [f"{cond}_rep{r}" for cond in CONDITIONS for r in range(1, N_REPS + 1)]

    # baseline mean expression per gene (log-normal)
    base_mu = RNG.lognormal(mean=4.0, sigma=1.5, size=N_GENES)

    # build count matrix
    counts = {}
    de_up = set(range(0, N_DE // 2))
    de_down = set(range(N_DE // 2, N_DE))
    de_genes = de_up | de_down

    for cond in CONDITIONS:
        mu = base_mu.copy()
        if cond == "treatment":
            mu[list(de_up)] *= DE_FOLD
            mu[list(de_down)] /= DE_FOLD
        for r in range(N_REPS):
            sample = f"{cond}_rep{r + 1}"
            # add biological variability then Poisson-sample
            lam = mu * RNG.gamma(10, 0.1, N_GENES)  # gene-wise noise
            counts[sample] = RNG.poisson(lam)

    mat = pd.DataFrame(counts, index=genes)
    mat.index.name = "gene"
    mat.to_csv(OUT / "counts.tsv", sep="\t")

    # truth table
    truth = pd.DataFrame({
        "gene": genes,
        "is_de": [i in de_genes for i in range(N_GENES)],
        "direction": ["up" if i in de_up else "down" if i in de_down else "none"
                       for i in range(N_GENES)],
    })
    truth.to_csv(OUT / "truth.tsv", sep="\t", index=False)
    print(f"wrote {mat.shape[0]} genes x {mat.shape[1]} samples -> {OUT / 'counts.tsv'}")
    print(f"planted DE genes: {N_DE} ({N_DE // 2} up, {N_DE // 2} down, fold={DE_FOLD})")
