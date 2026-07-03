#!/usr/bin/env python3
"""Synthesize demultiplexing reads with index hopping for the demux QC demo.

Models a patterned flow-cell run with 8 samples, each identified by a unique
(i7, i5) dual-index pair.  A small fraction of reads (~2%) undergo index
hopping (one index swapped to another sample's), and ~1% are undetermined
(neither index matches any known combo).

All data is synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(1)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

N_SAMPLES = 8
TOTAL_READS = 50_000
HOP_RATE = 0.02          # fraction of reads with one swapped index
UNDET_RATE = 0.01        # fraction with unknown indices

# generate unique i7 / i5 index sequences (8 bp each)
BASES = list("ACGT")
def _rand_index():
    return "".join(RNG.choice(BASES, 8))

i7_indices = [_rand_index() for _ in range(N_SAMPLES)]
i5_indices = [_rand_index() for _ in range(N_SAMPLES)]

if __name__ == "__main__":
    # per-sample read counts drawn from a Dirichlet so libraries aren't equal
    fracs = RNG.dirichlet(np.ones(N_SAMPLES) * 5)
    sample_counts = np.round(fracs * TOTAL_READS).astype(int)
    sample_counts[-1] = TOTAL_READS - sample_counts[:-1].sum()  # fix rounding

    rows = []
    for s in range(N_SAMPLES):
        n = sample_counts[s]
        for _ in range(n):
            r = RNG.random()
            if r < UNDET_RATE:
                # undetermined: random indices not in our known set
                i7, i5 = _rand_index(), _rand_index()
            elif r < UNDET_RATE + HOP_RATE:
                # index hop: keep one correct, swap the other
                if RNG.random() < 0.5:
                    i7 = i7_indices[s]
                    others = [j for j in range(N_SAMPLES) if j != s]
                    i5 = i5_indices[RNG.choice(others)]
                else:
                    i5 = i5_indices[s]
                    others = [j for j in range(N_SAMPLES) if j != s]
                    i7 = i7_indices[RNG.choice(others)]
            else:
                i7, i5 = i7_indices[s], i5_indices[s]
            rows.append((i7, i5))

    df = pd.DataFrame(rows, columns=["i7", "i5"])
    df = df.sample(frac=1, random_state=2).reset_index(drop=True)
    df.index.name = "read_id"

    # save the sample sheet alongside the reads
    ss = pd.DataFrame({"sample": [f"S{i+1}" for i in range(N_SAMPLES)],
                        "i7": i7_indices, "i5": i5_indices})
    ss.to_csv(OUT / "sample_sheet.csv", index=False)
    df.to_csv(OUT / "reads.csv")
    print(f"wrote {len(df)} reads for {N_SAMPLES} samples -> {OUT / 'reads.csv'}")
    print(f"wrote sample sheet -> {OUT / 'sample_sheet.csv'}")
