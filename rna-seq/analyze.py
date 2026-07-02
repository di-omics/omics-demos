#!/usr/bin/env python3
"""Bulk RNA-seq mini-pipeline: CPM normalize, PCA, per-gene t-test with BH."""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"


def cpm(counts):
    """Counts per million."""
    lib_size = counts.sum(axis=0)
    return counts / lib_size * 1e6


def bh_correct(pvals):
    """Benjamini-Hochberg FDR correction."""
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = np.empty(n)
    ranked[order] = np.arange(1, n + 1)
    return np.minimum(1.0, pvals * n / ranked)


def ttest_per_gene(mat, group_a_cols, group_b_cols):
    """Welch's t-test per gene, returns log2FC and p-values."""
    a = mat[group_a_cols].values.astype(float)
    b = mat[group_b_cols].values.astype(float)
    mean_a, mean_b = a.mean(axis=1), b.mean(axis=1)
    var_a, var_b = a.var(axis=1, ddof=1), b.var(axis=1, ddof=1)
    na, nb = a.shape[1], b.shape[1]
    se = np.sqrt(var_a / na + var_b / nb)
    se[se == 0] = np.inf
    t_stat = (mean_b - mean_a) / se
    # degrees of freedom (Welch-Satterthwaite)
    num = (var_a / na + var_b / nb) ** 2
    denom = (var_a / na) ** 2 / (na - 1) + (var_b / nb) ** 2 / (nb - 1)
    denom[denom == 0] = 1
    df = num / denom
    # two-tailed p-value from t-distribution (use normal approx for simplicity)
    from math import erfc, sqrt
    pvals = np.array([erfc(abs(t) / sqrt(2)) for t in t_stat])
    log2fc = np.log2((mean_b + 1) / (mean_a + 1))
    return log2fc, t_stat, pvals


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t", index_col=0)
    truth = pd.read_csv(DATA / "truth.tsv", sep="\t")

    # CPM + log2 transform
    normed = np.log2(cpm(counts) + 1)

    # PCA (manual, no sklearn needed)
    X = normed.values.T  # samples x genes
    X_centered = X - X.mean(axis=0)
    cov = np.cov(X_centered, rowvar=True)  # samples x samples
    eigvals, eigvecs = np.linalg.eigh(cov)
    idx = np.argsort(eigvals)[::-1]
    eigvals, eigvecs = eigvals[idx], eigvecs[:, idx]
    pcs = eigvecs[:, :2]  # project into PC space (already in sample space)
    var_explained = eigvals[:2] / eigvals.sum() * 100

    pca_df = pd.DataFrame({
        "sample": counts.columns,
        "PC1": pcs[:, 0], "PC2": pcs[:, 1],
        "condition": ["control" if "control" in s else "treatment" for s in counts.columns],
    })
    pca_df.to_csv(DATA / "pca.tsv", sep="\t", index=False)

    # DE testing
    ctrl = [c for c in counts.columns if "control" in c]
    treat = [c for c in counts.columns if "treatment" in c]
    log2fc, t_stat, pvals = ttest_per_gene(normed, ctrl, treat)
    padj = bh_correct(pvals)

    de_df = pd.DataFrame({
        "gene": counts.index,
        "log2FC": log2fc, "t_stat": t_stat,
        "pvalue": pvals, "padj": padj,
    })
    de_df.to_csv(DATA / "de_results.tsv", sep="\t", index=False)

    # report recovery
    sig = de_df[de_df.padj < 0.05]
    planted = set(truth[truth.is_de].gene)
    recovered = set(sig.gene) & planted
    fp = set(sig.gene) - planted

    print("\n=== Bulk RNA-seq DE analysis (synthetic) ===")
    print(f"PC1 explains {var_explained[0]:.1f}%, PC2 explains {var_explained[1]:.1f}%")
    print(f"Significant genes (padj < 0.05): {len(sig)}")
    print(f"Planted DE genes: {len(planted)}")
    print(f"Recovered: {len(recovered)} / {len(planted)} ({100*len(recovered)/len(planted):.0f}%)")
    print(f"False positives: {len(fp)}")
    print(f"\nwrote {DATA / 'de_results.tsv'}")


if __name__ == "__main__":
    main()
