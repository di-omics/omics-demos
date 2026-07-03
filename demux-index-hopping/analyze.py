#!/usr/bin/env python3
"""Demux QC: assign reads by index pair, build hopping matrix, quantify rate.

Reads are assigned to a sample when both i7 and i5 match a known pair.
Off-diagonal entries in the i7 x i5 count matrix reveal index hopping.
"""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"


def main():
    reads = pd.read_csv(DATA / "reads.csv")
    ss = pd.read_csv(DATA / "sample_sheet.csv")

    # build lookup: (i7, i5) -> sample name
    pair_to_sample = {(r.i7, r.i5): r.sample for _, r in ss.iterrows()}
    i7_list = ss.i7.tolist()
    i5_list = ss.i5.tolist()

    # assign each read
    reads["assigned"] = reads.apply(
        lambda r: pair_to_sample.get((r.i7, r.i5), "undetermined"), axis=1
    )

    # i7 x i5 count matrix (only known indices)
    i7_idx = {v: i for i, v in enumerate(i7_list)}
    i5_idx = {v: i for i, v in enumerate(i5_list)}
    n = len(i7_list)
    matrix = np.zeros((n, n), dtype=int)
    for _, r in reads.iterrows():
        ri = i7_idx.get(r.i7)
        ci = i5_idx.get(r.i5)
        if ri is not None and ci is not None:
            matrix[ri, ci] += 1

    sample_names = ss["sample"].tolist()
    mat_df = pd.DataFrame(matrix, index=sample_names, columns=sample_names)
    mat_df.index.name = "i7_sample"
    mat_df.columns.name = "i5_sample"

    diag_total = int(np.trace(matrix))
    off_diag_total = int(matrix.sum() - diag_total)
    classified = int(matrix.sum())
    undetermined = int((reads["assigned"] == "undetermined").sum())
    hop_rate = off_diag_total / classified if classified else 0

    # per-sample summary
    per_sample = pd.DataFrame({
        "sample": sample_names,
        "assigned_reads": [int(matrix[i, i]) for i in range(n)],
    })

    # flag unexpected combos
    flagged = []
    for i in range(n):
        for j in range(n):
            if i != j and matrix[i, j] > 0:
                flagged.append({
                    "i7_from": sample_names[i],
                    "i5_from": sample_names[j],
                    "count": int(matrix[i, j]),
                })
    flag_df = pd.DataFrame(flagged)

    print("\n=== Demux index-hopping QC (synthetic) ===")
    print(f"Total reads:           {len(reads)}")
    print(f"Classified (known i7xi5): {classified}")
    print(f"  On-diagonal (correct):  {diag_total}")
    print(f"  Off-diagonal (hopped):  {off_diag_total}")
    print(f"Undetermined:          {undetermined}")
    print(f"Estimated hopping rate:{hop_rate:8.2%}")
    print("\nPer-sample assigned reads:")
    print(per_sample.to_string(index=False))
    print("\ni7 x i5 count matrix:")
    print(mat_df.to_string())

    mat_df.to_csv(DATA / "hopping_matrix.csv")
    per_sample.to_csv(DATA / "per_sample.csv", index=False)
    if not flag_df.empty:
        flag_df.to_csv(DATA / "flagged_combos.csv", index=False)
    print(f"\nwrote {DATA / 'hopping_matrix.csv'} and {DATA / 'per_sample.csv'}")


if __name__ == "__main__":
    main()
