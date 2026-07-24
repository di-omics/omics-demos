#!/usr/bin/env python3
"""scRNA-seq UMI QC: three methods plus library-complexity rarefaction.

Dedup methods:
  1. Naive - unique (gene, UMI) pairs
  2. Directional - 1-mismatch merge, high-count parent absorbs low-count child
  3. Adjacency - 1-mismatch connected-component merge (UMI-tools-style)

The three diverge: naive overcounts (seq errors mint false UMIs), directional
is conservative, adjacency merges more aggressively.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations

DATA = Path(__file__).parent / "data"


def hamming1(a, b):
    return sum(x != y for x, y in zip(a, b)) == 1


def collapse_directional(umi_counts):
    """Directional 1-mismatch collapse. Returns #molecules."""
    umis = sorted(umi_counts, key=lambda u: -umi_counts[u])
    parent = {u: u for u in umis}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]; u = parent[u]
        return u

    for a, b in combinations(umis, 2):
        if hamming1(a, b):
            hi, lo = (a, b) if umi_counts[a] >= umi_counts[b] else (b, a)
            if umi_counts[hi] >= 2 * umi_counts[lo] - 1:
                parent[find(lo)] = find(hi)
    return len({find(u) for u in umis})


def collapse_adjacency(umi_counts):
    """Adjacency 1-mismatch collapse (UMI-tools style). Returns #molecules.

    Any two UMIs within Hamming distance 1 are merged into the same component,
    regardless of count ratio. More aggressive than directional.
    """
    umis = list(umi_counts.keys())
    parent = {u: u for u in umis}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]; u = parent[u]
        return u

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in combinations(umis, 2):
        if hamming1(a, b):
            union(a, b)
    return len({find(u) for u in umis})


def main():
    reads = pd.read_csv(DATA / "reads.tsv", sep="\t")
    total_reads = len(reads)

    raw = reads.groupby("gene").size().rename("raw_reads")
    naive = reads.drop_duplicates(["gene", "umi"]).groupby("gene").size().rename("naive_umi")

    dir_coll, adj_coll = {}, {}
    for gene, g in reads.groupby("gene"):
        uc = g.umi.value_counts().to_dict()
        dir_coll[gene] = collapse_directional(uc)
        adj_coll[gene] = collapse_adjacency(uc)
    dir_s = pd.Series(dir_coll, name="directional_umi")
    adj_s = pd.Series(adj_coll, name="adjacency_umi")

    tab = pd.concat([raw, naive, dir_s, adj_s], axis=1).fillna(0).astype(int)
    tab = tab.sort_values("raw_reads", ascending=False)
    dup_rate = 1 - tab.directional_umi.sum() / total_reads

    print("\n=== scRNA-seq UMI deduplication QC (synthetic) ===")
    print(f"Total reads:              {total_reads}")
    print(f"Naive-unique molecules:   {tab.naive_umi.sum()}")
    print(f"Directional (1-mm):       {tab.directional_umi.sum()}")
    print(f"Adjacency (1-mm):         {tab.adjacency_umi.sum()}")
    print(f"Duplication rate (dir):    {dup_rate:.1%}")
    print("\nTop genes:")
    print(tab.head(8).to_string())

    # library-complexity rarefaction: subsample reads, count molecules by each method
    order = reads.sample(frac=1, random_state=3).reset_index(drop=True)
    depths = np.linspace(0.05, 1.0, 20)
    sat = []
    for d in depths:
        n = int(d * total_reads)
        sub = order.iloc[:n]
        n_naive = sub.drop_duplicates(["gene", "umi"]).shape[0]
        n_dir, n_adj = 0, 0
        for _, g in sub.groupby("gene"):
            uc = g.umi.value_counts().to_dict()
            n_dir += collapse_directional(uc)
            n_adj += collapse_adjacency(uc)
        sat.append((n, n_naive, n_dir, n_adj))
    sat_df = pd.DataFrame(sat, columns=["reads", "naive", "directional", "adjacency"])
    sat_df.to_csv(DATA / "saturation.tsv", sep="\t", index=False)
    tab.to_csv(DATA / "counts.tsv", sep="\t")
    print(f"\nwrote {DATA / 'counts.tsv'} and {DATA / 'saturation.tsv'}")


if __name__ == "__main__":
    main()
