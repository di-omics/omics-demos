#!/usr/bin/env python3
"""Compute EM-seq validation QC from synthetic methylation calls.

Metrics:
  - Conversion efficiency   = 1 - mCpG(lambda)     [target > 99.5%]
  - CpG protection          = mCpG(pUC19)          [target > ~95%]
  - Global CpG methylation  per condition          [human ~70-80%]
  - Mean coverage           per condition
  - 10 ng vs 0.1 ng comparison
"""
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"


def meth_frac(df):
    m, u = df.meth_reads.sum(), df.unmeth_reads.sum()
    return m / (m + u)


def load(name):
    return pd.read_csv(DATA / name, sep="\t")


def main():
    conv_eff = 1 - meth_frac(load("control_lambda.tsv"))
    protection = meth_frac(load("control_puc19.tsv"))

    rows = []
    for cond, f in [("10 ng", "cpg_10ng.tsv"), ("0.1 ng", "cpg_0p1ng.tsv")]:
        d = load(f)
        cov = d.meth_reads + d.unmeth_reads
        rows.append({
            "condition": cond,
            "cpgs_covered": len(d),
            "mean_coverage": round(cov.mean(), 2),
            "global_mCpG": round(meth_frac(d), 4),
        })
    summary = pd.DataFrame(rows)

    print("\n=== EM-seq validation QC (synthetic) ===")
    print(f"Conversion efficiency (lambda): {conv_eff:.2%}   [target > 99.5%]")
    print(f"CpG protection (pUC19):         {protection:.2%}   [target > 95%]")
    print("\nPer-condition:")
    print(summary.to_string(index=False))

    out = DATA / "qc_summary.tsv"
    summary.assign(
        conversion_efficiency=round(conv_eff, 5),
        puc19_protection=round(protection, 5),
    ).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
