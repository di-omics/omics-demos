# emseq-methylation-demo

Toy EM-seq methylation QC pipeline. Takes per-CpG methylation calls plus spike-in
conversion controls and reports the metrics you actually gate an EM-seq run on -
conversion efficiency, CpG protection, coverage, and global methylation - across a
standard (10 ng) and an ultra-low (0.1 ng) input.

All data is synthetic. This is a self-contained demo of the *analysis* pattern, not
a production pipeline.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make emseq
```

Generates synthetic data, runs QC, and writes a plot to `assets/`.

## What it computes

- **Conversion efficiency** from an unmethylated lambda spike-in - `1 − mCpG(λ)`. Target > 99.5%.
- **CpG protection** from a CpG-methylated pUC19 spike-in - `mCpG(pUC19)`. Target > 95%.
- **Global CpG methylation** and **mean coverage**, per input condition.
- A **10 ng vs 0.1 ng** comparison - ultra-low input drops coverage and adds noise
  while conversion holds.

## Example output

```
=== EM-seq validation QC (synthetic) ===
Conversion efficiency (lambda): 99.71%   [target > 99.5%]
CpG protection (pUC19):         97.22%   [target > 95%]

Per-condition:
condition  cpgs_covered  mean_coverage  global_mCpG
    10 ng         20000          29.98        0.782
   0.1 ng         19947           5.99        0.770
```

![QC plots](assets/qc_plots.png)

## Layout

```
generate_data.py   synthesize CpG calls + spike-in controls
analyze.py         conversion efficiency, coverage, global methylation
plots.py           per-CpG methylation histogram + control bars
```

Run from the repo root with `make emseq`.

## License

MIT
