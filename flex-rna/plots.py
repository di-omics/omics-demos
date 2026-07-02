#!/usr/bin/env python3
"""Plots for the Flex demo (writes assets/flex_qc.png)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)


def main():
    df = pd.read_csv(DATA / "cell_qc.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))

    # barcode-rank (knee) plot
    ranked = df.total_umi.sort_values(ascending=False).values
    ax[0].loglog(range(1, len(ranked) + 1), ranked, color="#3b6ea5")
    thr = df[df.called_cell].total_umi.min()
    ax[0].axhline(thr, ls="--", color="#a53b6e", lw=1, label="cell threshold")
    ax[0].set(title="Barcode rank (knee)", xlabel="barcode rank", ylabel="total UMI")
    ax[0].legend(frameon=False)

    # probes/cell per sample
    called = df[df.called_cell]
    samples = sorted(called.sample_bc.unique())
    data = [called[called.sample_bc == s].probes_detected.values for s in samples]
    ax[1].boxplot(data, showfliers=False)
    ax[1].set_xticks(range(1, len(samples) + 1)); ax[1].set_xticklabels(samples)
    ax[1].set(title="Probes detected / cell", ylabel="probes", xlabel="sample barcode")
    fig.tight_layout(); out = ASSETS / "flex_qc.png"; fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
