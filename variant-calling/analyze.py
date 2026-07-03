#!/usr/bin/env python3
"""Variant calling from piled-up reads: call SNVs and score against truth.

Piles up reads per genomic position, computes per-base allele frequencies,
calls variants above a minimum AF threshold with sufficient depth, then
evaluates precision / recall against the planted truth set. Reports per-AF-bin
sensitivity so you can see recall drop at low allele fractions.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from collections import Counter

DATA = Path(__file__).parent / "data"

AF_THRESHOLD = 0.02      # minimum alt allele fraction to call a variant
MIN_DEPTH = 10            # minimum read depth to consider a position
AF_BIN_EDGES = [0.0, 0.05, 0.10, 0.20, 0.35, 0.60, 1.01]


def pileup(reads_df, ref_len):
    """Return per-position base counts as a list of Counters."""
    counts = [Counter() for _ in range(ref_len)]
    for _, row in reads_df.iterrows():
        start = row.start_pos
        for offset, base in enumerate(row.sequence):
            counts[start + offset][base] += 1
    return counts


def call_variants(counts, ref_bases):
    """Call variants where alt AF > threshold and depth > min_depth."""
    calls = []
    for pos, ctr in enumerate(counts):
        depth = sum(ctr.values())
        if depth < MIN_DEPTH:
            continue
        ref_base = ref_bases[pos]
        for base, count in ctr.items():
            if base == ref_base:
                continue
            af = count / depth
            if af > AF_THRESHOLD:
                calls.append((pos, ref_base, base, round(af, 4), depth))
    return pd.DataFrame(calls, columns=["pos", "ref", "alt", "obs_af", "depth"])


def score(calls_df, truth_df):
    """Score calls against truth; return TP, FP, FN sets and merged DataFrame."""
    truth_set = set(zip(truth_df.pos, truth_df.alt))
    call_set = set(zip(calls_df.pos, calls_df.alt))

    tp = truth_set & call_set
    fp = call_set - truth_set
    fn = truth_set - call_set

    precision = len(tp) / len(call_set) if call_set else 0
    recall = len(tp) / len(truth_set) if truth_set else 0

    # merge for per-variant detail
    merged = truth_df.merge(calls_df[["pos", "alt", "obs_af", "depth"]],
                            on=["pos", "alt"], how="left")
    merged["called"] = merged.obs_af.notna()
    return precision, recall, len(tp), len(fp), len(fn), merged


def sensitivity_by_bin(merged, edges):
    """Group true variants by AF bin and compute per-bin sensitivity."""
    merged = merged.copy()
    labels = [f"{edges[i]:.2f}-{edges[i+1]:.2f}" for i in range(len(edges) - 1)]
    merged["af_bin"] = pd.cut(merged.true_af, bins=edges, labels=labels,
                              include_lowest=True)
    grouped = merged.groupby("af_bin", observed=False).agg(
        total=("called", "size"),
        detected=("called", "sum"),
    )
    grouped["sensitivity"] = (grouped.detected / grouped.total).fillna(0).round(3)
    return grouped.reset_index()


def main():
    reads = pd.read_csv(DATA / "reads.tsv", sep="\t")
    truth = pd.read_csv(DATA / "truth.tsv", sep="\t")
    ref = pd.read_csv(DATA / "reference.tsv", sep="\t")
    ref_len = len(ref)

    print("piling up reads …")
    counts = pileup(reads, ref_len)

    ref_bases = ref.base.tolist()

    print("calling variants …")
    calls = call_variants(counts, ref_bases)

    precision, recall, n_tp, n_fp, n_fn, merged = score(calls, truth)

    print("\n=== Variant calling QC (synthetic) ===")
    print(f"True variants:   {len(truth)}")
    print(f"Called variants: {len(calls)}")
    print(f"TP / FP / FN:    {n_tp} / {n_fp} / {n_fn}")
    print(f"Precision:       {precision:.3f}")
    print(f"Recall:          {recall:.3f}")

    bin_df = sensitivity_by_bin(merged, AF_BIN_EDGES)
    print("\nSensitivity by allele-fraction bin:")
    print(bin_df.to_string(index=False))

    # save outputs
    calls.to_csv(DATA / "calls.tsv", sep="\t", index=False)
    merged.to_csv(DATA / "truth_vs_calls.tsv", sep="\t", index=False)
    bin_df.to_csv(DATA / "sensitivity_by_bin.tsv", sep="\t", index=False)
    print(f"\nwrote {DATA / 'calls.tsv'}, {DATA / 'truth_vs_calls.tsv'}, "
          f"{DATA / 'sensitivity_by_bin.tsv'}")


if __name__ == "__main__":
    main()
