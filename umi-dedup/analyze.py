#!/usr/bin/env python3
"""UMI dedup QC: raw reads vs collapsed molecules, duplication rate, saturation.

Collapsing uses a simple 1-mismatch directional merge within each gene: UMIs
within Hamming distance 1 are folded into the more abundant UMI. This recovers
true molecule counts from PCR-amplified, error-containing reads.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from itertools import combinations

DATA = Path(__file__).parent / "data"


def hamming1(a, b):
    return sum(x != y for x, y in zip(a, b)) == 1


def collapse_gene(umi_counts):
    """Directional 1-mismatch collapse. umi_counts: dict umi->count. Returns #molecules."""
    umis = sorted(umi_counts, key=lambda u: -umi_counts[u])
    parent = {u: u for u in umis}

    def find(u):
        while parent[u] != u:
            parent[u] = parent[parent[u]]; u = parent[u]
        return u

    for a, b in combinations(umis, 2):
        if hamming1(a, b):
            hi, lo = (a, b) if umi_counts[a] >= umi_counts[b] else (b, a)
            # only merge low-count child into higher-count parent (directional)
            if umi_counts[hi] >= 2 * umi_counts[lo] - 1:
                parent[find(lo)] = find(hi)
    return len({find(u) for u in umis})


def main():
    reads = pd.read_csv(DATA / "reads.tsv", sep="\t")
    total_reads = len(reads)

    raw = reads.groupby("gene").size().rename("raw_reads")
    # naive dedup = unique (gene, umi)
    naive = reads.drop_duplicates(["gene", "umi"]).groupby("gene").size().rename("naive_umi")
    # 1-mismatch collapsed
    coll = {}
    for gene, g in reads.groupby("gene"):
        coll[gene] = collapse_gene(g.umi.value_counts().to_dict())
    coll = pd.Series(coll, name="collapsed_umi")

    tab = pd.concat([raw, naive, coll], axis=1).fillna(0).astype(int).sort_values("raw_reads", ascending=False)
    dup_rate = 1 - tab.collapsed_umi.sum() / total_reads

    print("\n=== UMI dedup QC (synthetic) ===")
    print(f"Total reads:           {total_reads}")
    print(f"Naive-unique molecules:{tab.naive_umi.sum()}")
    print(f"Collapsed molecules:   {tab.collapsed_umi.sum()}   (1-mismatch merge)")
    print(f"Duplication rate:      {dup_rate:.1%}")
    print("\nTop genes:")
    print(tab.head(8).to_string())

    # saturation: unique molecules vs reads subsampled
    order = reads.sample(frac=1, random_state=3).reset_index(drop=True)
    depths = np.linspace(0.1, 1.0, 10)
    sat = []
    for d in depths:
        sub = order.iloc[: int(d * total_reads)]
        mols = sub.drop_duplicates(["gene", "umi"]).shape[0]
        sat.append((int(d * total_reads), mols))
    pd.DataFrame(sat, columns=["reads", "molecules"]).to_csv(DATA / "saturation.tsv", sep="\t", index=False)
    tab.to_csv(DATA / "counts.tsv", sep="\t")
    print(f"\nwrote {DATA/'counts.tsv'} and {DATA/'saturation.tsv'}")


if __name__ == "__main__":
    main()
