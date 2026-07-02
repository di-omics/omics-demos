# rna-seq

A small bulk RNA-seq pipeline: synthetic count matrix in, PCA and a differential
expression call out, scored against a known planted truth so you can see how much
signal the test actually recovers.

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
