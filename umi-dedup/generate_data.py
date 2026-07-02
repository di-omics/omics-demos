#!/usr/bin/env python3
"""Synthesize UMI-tagged reads for the dedup demo.

Models a FLASH-seq-style library: a set of true mRNA molecules, each tagged with
a random UMI, then PCR-amplified into many reads. A small per-base error rate
introduces near-duplicate UMIs (the thing UMI collapsing has to handle).

All data is synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(1)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

N_GENES = 40
N_MOLECULES = 3000          # true captured molecules
UMI_LEN = 10
PCR_MEAN_DUP = 6            # mean reads per molecule after amplification
UMI_ERR = 0.004            # per-base sequencing error

BASES = np.array(list("ACGT"))


def random_umi():
    return "".join(RNG.choice(BASES, UMI_LEN))


def mutate(umi):
    s = list(umi)
    for i in range(len(s)):
        if RNG.random() < UMI_ERR:
            s[i] = RNG.choice(BASES)
    return "".join(s)


if __name__ == "__main__":
    # true molecules: gene assignment skewed so a few genes are highly expressed
    gene_p = RNG.dirichlet(np.ones(N_GENES) * 0.3)
    genes = RNG.choice(N_GENES, N_MOLECULES, p=gene_p)
    umis = [random_umi() for _ in range(N_MOLECULES)]

    rows = []
    rid = 0
    for g, u in zip(genes, umis):
        ndup = 1 + RNG.poisson(PCR_MEAN_DUP - 1)
        for _ in range(ndup):
            rows.append((f"r{rid}", f"gene_{g:02d}", mutate(u)))
            rid += 1
    df = pd.DataFrame(rows, columns=["read_id", "gene", "umi"])
    df = df.sample(frac=1, random_state=2).reset_index(drop=True)  # shuffle
    df.to_csv(OUT / "reads.tsv", sep="\t", index=False)
    print(f"wrote {len(df)} reads across {N_MOLECULES} true molecules -> {OUT/'reads.tsv'}")
