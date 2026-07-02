# rna-seq

A small bulk RNA-seq pipeline: synthetic count matrix in, PCA and a differential
expression call out, scored against a known planted truth so you can see how much
signal the test actually recovers.

All data is synthetic - no real samples.

## Run

```bash
pip install -r requirements.txt   # includes scipy
make rna
```

## What it does

1. **Simulate** - 2,000 genes across two conditions (5 replicates each). Gene
   baselines are log-normal; a planted set of 120 genes gets a real log2 fold-change
   in the treated arm. Counts are negative-binomial with per-sample library-size
   variation, so normalization matters.
2. **Normalize** - counts to CPM, then log2(CPM+1).
3. **PCA** - on the top-variable genes (numpy SVD), to check the samples separate by
   condition rather than by depth or noise.
4. **Test** - per-gene Welch t-test, Benjamini-Hochberg FDR, and call DE at
   FDR < 0.05 with |log2FC| >= 1.
5. **Score** - recall, precision, and sign agreement against the planted truth.

## Example output

```
=== RNA-seq differential expression (synthetic) ===
2000 genes x 10 samples  |  1970 genes pass expression filter
PCA: PC1 42.5% var, PC2 9.9% var
called DE: 78 genes at FDR<0.05, |log2FC|>=1.0
vs 120 planted truth -> recall 57.5%, precision 88.5%, sign agreement 100.0%
```

PC1 captures the condition split, and the test recovers most of the planted DE at
high precision - the misses are the small-effect genes that n=5 can't separate from
noise, which is the honest outcome.

![PCA and volcano](assets/rnaseq_qc.png)

## How it works

Per-gene Welch t-test with Benjamini-Hochberg FDR correction (from `analyze.py`):

```python
def bh(p):
    p = np.asarray(p); n = len(p); order = np.argsort(p)
    q = np.empty(n); q[order] = (p[order] * n / (np.arange(n) + 1))
    ranked = q[order]
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    q[order] = np.clip(ranked, 0, 1)
    return q

# per-gene Welch t-test, ctrl vs treat
ctrl = meta.query("condition=='ctrl'")["sample"].values
treat = meta.query("condition=='treat'")["sample"].values
a, b = logf[ctrl].values, logf[treat].values
t, p = stats.ttest_ind(b, a, axis=1, equal_var=False)
lfc = b.mean(1) - a.mean(1)
res["padj"] = bh(res.pval.values)
res["called"] = (res.padj < FDR) & (res.log2fc.abs() >= LFC)
```

Genes are called DE at FDR < 0.05 with |log2FC| >= 1, then scored against the planted truth.

## Files

```
generate_data.py   simulate the count matrix + planted truth
analyze.py         CPM/log normalize, PCA, Welch t-test + BH FDR, score vs truth
plots.py           sample PCA + volcano
```

## Outputs

```
data/de_results.tsv   per-gene log2FC, p, adjusted p, called flag
data/pca.tsv          per-sample PC1/PC2 + condition
assets/rnaseq_qc.png  PCA scatter + volcano
```
