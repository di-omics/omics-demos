#!/usr/bin/env python3
"""Compute methylation-sequencing QC from synthetic methylation calls.

Metrics:
  - Conversion efficiency from an unmethylated control
  - CpG protection from a methylated control
  - Global methylation per context (CpG, CHG, CHH) per condition
  - Mean coverage per condition
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
    conv_eff = 1 - meth_frac(load("control_unmethylated.tsv"))
    protection = meth_frac(load("control_methylated.tsv"))

    rows = []
    for cond, f in [("10 ng", "meth_10ng.tsv"), ("0.1 ng", "meth_0p1ng.tsv")]:
        d = load(f)
        cov = d.meth_reads + d.unmeth_reads
        cpg = d[d.context == "CpG"]
        chg = d[d.context == "CHG"]
        chh = d[d.context == "CHH"]
        rows.append({
            "condition": cond,
            "sites_total": len(d),
            "mean_coverage": round(cov.mean(), 2),
            "mCpG": round(meth_frac(cpg), 4),
            "mCHG": round(meth_frac(chg), 4),
            "mCHH": round(meth_frac(chh), 4),
        })
    summary = pd.DataFrame(rows)

    # --- validation against planted truth ---
    TRUE_MCPG = 0.78
    TRUE_UNMETHYLATED_CONTROL = 0.003
    TRUE_METHYLATED_CONTROL = 0.972

    print("\n=== Methylation-sequencing QC (synthetic) ===")
    print(f"Conversion efficiency: {conv_eff:.2%}   [target > 99.5%]")
    print(f"CpG protection:         {protection:.2%}   [target > 95%]")
    print("\nPer-condition methylation by context:")
    print(summary.to_string(index=False))
    print(f"\nRecovery vs planted truth:")
    for _, r in summary.iterrows():
        print(f"  {r['condition']}: mCpG {r['mCpG']:.3f} (truth {TRUE_MCPG})")
    print(f"  unmethylated-control rate {1 - conv_eff:.4f} "
          f"(truth {TRUE_UNMETHYLATED_CONTROL})")
    print(f"  methylated-control rate   {protection:.4f} "
          f"(truth {TRUE_METHYLATED_CONTROL})")

    out = DATA / "qc_summary.tsv"
    summary.assign(
        conversion_efficiency=round(conv_eff, 5),
        methylated_control_protection=round(protection, 5),
    ).to_csv(out, sep="\t", index=False)
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
