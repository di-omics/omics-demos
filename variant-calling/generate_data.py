#!/usr/bin/env python3
"""Synthesize a reference genome, plant SNVs, and generate piled-up reads.

Creates a small reference (~10 kbp), plants ~30 SNVs at known positions with
varying allele fractions (0.05-1.0), then tiles ~200x coverage of 150 bp reads.
Each read flips variant bases with probability = true AF and adds a low
sequencing-error background at non-variant positions.

All data is synthetic.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(1)
OUT = Path(__file__).parent / "data"; OUT.mkdir(exist_ok=True)

REF_LEN = 10_000
N_VARIANTS = 30
READ_LEN = 150
TARGET_DEPTH = 200
SEQ_ERR = 0.005          # 0.5 % per-base sequencing error at non-variant sites
BASES = list("ACGT")

if __name__ == "__main__":
    # --- reference sequence ---
    ref = RNG.choice(BASES, REF_LEN)

    # --- plant variants ---
    var_positions = np.sort(RNG.choice(REF_LEN, N_VARIANTS, replace=False))
    # spread AFs from very low to 1.0 so sensitivity drops at the low end
    true_af = np.concatenate([
        RNG.uniform(0.01, 0.04, 6),     # 6 very-low-AF variants
        RNG.uniform(0.04, 0.10, 5),      # 5 low-AF variants
        RNG.uniform(0.10, 0.35, 7),      # 7 moderate
        RNG.uniform(0.35, 1.0, N_VARIANTS - 18),  # rest high
    ])
    RNG.shuffle(true_af)
    alt_alleles = []
    for pos in var_positions:
        others = [b for b in BASES if b != ref[pos]]
        alt_alleles.append(RNG.choice(others))

    truth_rows = []
    for pos, alt, af in zip(var_positions, alt_alleles, true_af):
        truth_rows.append((int(pos), ref[pos], alt, round(float(af), 4)))
    truth_df = pd.DataFrame(truth_rows, columns=["pos", "ref", "alt", "true_af"])
    truth_df.to_csv(OUT / "truth.tsv", sep="\t", index=False)

    # variant lookup for fast read generation
    var_lookup = {int(pos): (alt, af) for pos, alt, af in
                  zip(var_positions, alt_alleles, true_af)}

    # --- reference.tsv ---
    ref_df = pd.DataFrame({"pos": range(REF_LEN), "base": ref})
    ref_df.to_csv(OUT / "reference.tsv", sep="\t", index=False)

    # --- generate reads ---
    n_reads = int(TARGET_DEPTH * REF_LEN / READ_LEN)
    starts = RNG.integers(0, REF_LEN - READ_LEN, n_reads)

    read_rows = []
    for rid, start in enumerate(starts):
        seq = list(ref[start:start + READ_LEN])
        for offset in range(READ_LEN):
            gpos = start + offset
            if gpos in var_lookup:
                alt, af = var_lookup[gpos]
                if RNG.random() < af:
                    seq[offset] = alt
            else:
                if RNG.random() < SEQ_ERR:
                    others = [b for b in BASES if b != seq[offset]]
                    seq[offset] = RNG.choice(others)
        read_rows.append((f"r{rid}", int(start), "".join(seq)))

    reads_df = pd.DataFrame(read_rows, columns=["read_id", "start_pos", "sequence"])
    reads_df.to_csv(OUT / "reads.tsv", sep="\t", index=False)

    print(f"wrote {len(ref_df)} ref bases   -> {OUT / 'reference.tsv'}")
    print(f"wrote {len(truth_df)} true SNVs  -> {OUT / 'truth.tsv'}")
    print(f"wrote {len(reads_df)} reads      -> {OUT / 'reads.tsv'}")
