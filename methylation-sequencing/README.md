# Methylation-sequencing QC

A synthetic methylation-sequencing QC pipeline. It evaluates conversion
efficiency, CpG protection, coverage, global methylation, and read-position bias
across standard and low-input conditions.

The demo uses functional unmethylated and methylated controls. All data is
generated locally.

## Run

```bash
pip install -r requirements.txt
make methylation
```

## Metrics

- Conversion efficiency from an unmethylated control
- CpG protection from a methylated control
- Global CpG, CHG, and CHH methylation
- Mean coverage by input condition
- M-bias across read positions

## Example output

```text
=== Methylation-sequencing QC (synthetic) ===
Conversion efficiency: 99.69%   [target > 99.5%]
CpG protection:         97.24%   [target > 95%]

condition  sites_total  mean_coverage  mCpG   mCHG   mCHH
    10 ng        75000          30.00 0.782 0.0032 0.0032
   0.1 ng        74799           6.02 0.775 0.0042 0.0041
```

Expected result: recovered methylation remains close to the planted truth while
the low-input condition loses coverage and gains noise.

![Methylation-sequencing QC](assets/qc_plots.png)

## Files

```text
generate_data.py   synthetic context calls, controls, and M-bias
analyze.py         conversion, protection, coverage, and methylation metrics
plots.py           methylation distributions, controls, contexts, and M-bias
```

All values are illustrative and are not a validated clinical threshold set.
