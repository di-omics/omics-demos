#!/usr/bin/env python3
"""Normalize, run PCA, and call differential expression on the synthetic RNA-seq
counts, then score the calls against the planted truth.

  CPM -> log2(CPM+1) -> PCA (numpy SVD) -> per-gene Welch t-test -> BH FDR.

Writes de_results.tsv and pca.tsv.
"""
import numpy as np, pandas as pd
from scipy import stats
from pathlib import Path

DATA = Path(__file__).parent / "data"
MIN_MEAN, FDR, LFC = 1.0, 0.05, 1.0


def bh(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    q = np.empty(n); q[order] = (p[order] * n / (np.arange(n) + 1))
    # enforce monotonicity from the largest p downward
    ranked = q[order]
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return q


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t", index_col=0)
    meta = pd.read_csv(DATA / "samples.tsv", sep="\t")
    truth = pd.read_csv(DATA / "truth_de.tsv", sep="\t")

    # CPM + log
    cpm = counts / counts.sum(axis=0) * 1e6
    logcpm = np.log2(cpm + 1)
    keep = counts.mean(axis=1) >= MIN_MEAN
    logf = logcpm[keep]

    # PCA on top-variable genes
    top = logf.var(axis=1).sort_values(ascending=False).head(500).index
    X = logf.loc[top].T.values
    X = X - X.mean(0)
    U, Sv, Vt = np.linalg.svd(X, full_matrices=False)
    scores = U * Sv
    var_exp = (Sv ** 2) / (Sv ** 2).sum()
    pca = pd.DataFrame({"sample": logf.columns, "condition": meta.set_index("sample").loc[logf.columns, "condition"].values,
                        "PC1": scores[:, 0], "PC2": scores[:, 1]})
    pca.to_csv(DATA / "pca.tsv", sep="\t", index=False)

    # per-gene Welch t-test, ctrl vs treat
    ctrl = meta.query("condition=='ctrl'")["sample"].values
    treat = meta.query("condition=='treat'")["sample"].values
    a, b = logf[ctrl].values, logf[treat].values
    t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)
    lfc = b.mean(1) - a.mean(1)
    res = pd.DataFrame({"gene": logf.index, "log2fc": lfc, "pval": p})
    res["padj"] = bh(res.pval.values)
    res["called"] = (res.padj < FDR) & (res.log2fc.abs() >= LFC)
    res = res.sort_values("padj")
    res.to_csv(DATA / "de_results.tsv", sep="\t", index=False)

    # score vs truth
    truth_set = set(truth.gene)
    called_set = set(res.query("called").gene)
    tp = len(called_set & truth_set)
    recall = tp / len(truth_set)
    prec = tp / max(len(called_set), 1)
    # directional agreement on recovered genes
    merged = res[res.gene.isin(truth_set)].merge(truth, on="gene")
    sign_ok = (np.sign(merged.log2fc) == np.sign(merged.log2fc_true)).mean()

    print("=== RNA-seq differential expression (synthetic) ===")
    print(f"{counts.shape[0]} genes x {counts.shape[1]} samples  |  {keep.sum()} genes pass expression filter")
    print(f"PCA: PC1 {100*var_exp[0]:.1f}% var, PC2 {100*var_exp[1]:.1f}% var")
    print(f"called DE: {len(called_set)} genes at FDR<{FDR}, |log2FC|>={LFC}")
    print(f"vs {len(truth_set)} planted truth -> recall {recall:.1%}, precision {prec:.1%}, "
          f"sign agreement {sign_ok:.1%}")
    print(f"\nwrote {DATA/'de_results.tsv'}, {DATA/'pca.tsv'}")


if __name__ == "__main__":
    main()
