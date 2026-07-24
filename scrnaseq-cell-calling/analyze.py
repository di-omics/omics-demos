#!/usr/bin/env python3
"""scRNA-seq QC: knee-point cell calling, ambient RNA, and sample summaries."""
import numpy as np
import pandas as pd
from pathlib import Path

DATA = Path(__file__).parent / "data"


def find_knee(ranked_umis):
    """Find the knee point via max-curvature on the log-log barcode-rank curve.

    Fits discrete curvature on the log-transformed rank vs UMI curve and returns
    the rank index with maximum curvature.
    """
    n = len(ranked_umis)
    x = np.log10(np.arange(1, n + 1))
    y = np.log10(ranked_umis.clip(min=1))
    # smooth with a moving average to avoid noise spikes
    win = max(5, n // 50)
    y_smooth = np.convolve(y, np.ones(win) / win, mode="same")
    # discrete curvature: |y''| / (1 + y'^2)^(3/2)
    dy = np.gradient(y_smooth, x)
    ddy = np.gradient(dy, x)
    curvature = np.abs(ddy) / (1 + dy ** 2) ** 1.5
    # ignore edges
    margin = max(10, n // 20)
    curvature[:margin] = 0
    curvature[-margin:] = 0
    knee_idx = np.argmax(curvature)
    return knee_idx, ranked_umis[knee_idx]


def estimate_ambient(probe_mat, is_empty):
    """Estimate ambient RNA profile from empty barcodes.

    Returns per-probe ambient fraction (sums to 1).
    """
    empty_counts = probe_mat[is_empty].sum(axis=0)
    total = empty_counts.sum()
    if total == 0:
        return np.zeros(probe_mat.shape[1])
    return (empty_counts / total).values


def main():
    counts = pd.read_csv(DATA / "counts.tsv", sep="\t")
    bc = pd.read_csv(DATA / "barcodes.tsv", sep="\t")
    probe_cols = [c for c in counts.columns if c.startswith("probe_")]
    m = counts[probe_cols]

    total_umi = m.sum(axis=1)
    n_probes = (m > 0).sum(axis=1)

    # knee-point cell calling via max curvature
    ranked = total_umi.sort_values(ascending=False).values
    knee_idx, knee_umi = find_knee(ranked)
    thresh = knee_umi
    is_cell = total_umi >= thresh

    df = pd.DataFrame({
        "barcode": counts.barcode, "sample_bc": bc.sample_bc, "truth": bc.truth,
        "total_umi": total_umi, "probes_detected": n_probes, "called_cell": is_cell,
    })

    # ambient RNA estimate from barcodes below threshold
    is_empty_called = ~is_cell
    ambient_profile = estimate_ambient(m, is_empty_called.values)
    ambient_frac = pd.DataFrame({
        "probe": probe_cols,
        "ambient_fraction": ambient_profile,
    })
    ambient_frac.to_csv(DATA / "ambient_profile.tsv", sep="\t", index=False)

    # ambient contamination estimate for cells
    empty_mean_umi = total_umi[is_empty_called].mean() if is_empty_called.any() else 0
    cell_mean_umi = total_umi[is_cell].mean() if is_cell.any() else 1
    ambient_pct = 100 * empty_mean_umi / cell_mean_umi

    print("\n=== scRNA-seq cell-calling QC (synthetic; max-curvature knee) ===")
    print(f"Knee point: rank {knee_idx + 1}, UMI threshold {thresh:.0f}")
    tp = ((df.truth == "cell") & df.called_cell).sum()
    fp = ((df.truth == "empty") & df.called_cell).sum()
    fn = ((df.truth == "cell") & ~df.called_cell).sum()
    print(f"Called cells: {df.called_cell.sum()}   (true cells kept {tp}, "
          f"empties admitted {fp}, cells lost {fn})")
    print(f"Ambient RNA estimate: {ambient_pct:.1f}% of cell UMI "
          f"(empty mean {empty_mean_umi:.0f}, cell mean {cell_mean_umi:.0f})")

    print("\nPer-sample (called cells):")
    called = df[df.called_cell]
    summ = called.groupby("sample_bc").agg(
        cells=("barcode", "size"),
        median_umi=("total_umi", "median"),
        median_probes=("probes_detected", "median"),
    )
    print(summ.to_string())

    df.to_csv(DATA / "cell_qc.tsv", sep="\t", index=False)
    print(f"\nwrote {DATA / 'cell_qc.tsv'}, {DATA / 'ambient_profile.tsv'}")


if __name__ == "__main__":
    main()
