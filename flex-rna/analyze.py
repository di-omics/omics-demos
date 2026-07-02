#!/usr/bin/env python3
"""Flex QC: cell calling (knee), per-sample demux summary, probes/cell."""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t")
    bc = pd.read_csv(DATA / "barcodes.tsv", sep="\t")
    probe_cols = [c for c in counts.columns if c.startswith("probe_")]
    m = counts[probe_cols]

    total_umi = m.sum(axis=1)
    n_probes = (m > 0).sum(axis=1)

    # simple knee-style cell call: UMI threshold at 10% of the 99th percentile
    thresh = np.percentile(total_umi, 99) * 0.10
    is_cell = total_umi >= thresh

    df = pd.DataFrame({
        "barcode": counts.barcode, "sample_bc": bc.sample_bc, "truth": bc.truth,
        "total_umi": total_umi, "probes_detected": n_probes, "called_cell": is_cell,
    })

    print("\n=== 10x Flex QC (synthetic; assumed Flex design) ===")
    print(f"UMI threshold for cell call: {thresh:.0f}")
    called = df[df.called_cell]
    # confusion vs truth
    tp = ((df.truth == "cell") & df.called_cell).sum()
    fp = ((df.truth == "empty") & df.called_cell).sum()
    fn = ((df.truth == "cell") & ~df.called_cell).sum()
    print(f"Called cells: {df.called_cell.sum()}   (true cells kept {tp}, empties admitted {fp}, cells lost {fn})")

    print("\nPer-sample (called cells):")
    summ = called.groupby("sample_bc").agg(
        cells=("barcode", "size"),
        median_umi=("total_umi", "median"),
        median_probes=("probes_detected", "median"),
    )
    print(summ.to_string())

    df.to_csv(DATA / "cell_qc.tsv", sep="\t", index=False)
    print(f"\nwrote {DATA/'cell_qc.tsv'}")


if __name__ == "__main__":
    main()
