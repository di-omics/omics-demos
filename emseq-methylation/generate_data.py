#!/usr/bin/env python3
"""Generate synthetic EM-seq methylation data for the demo.

Produces per-CpG methylation calls for two input conditions (10 ng and 0.1 ng)
plus unmethylated (lambda) and CpG-methylated (pUC19) spike-in controls, which
are the standard controls used to measure conversion efficiency in EM-seq.

All data is synthetic. No real sequencing data is used.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(0)
OUT = Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)

N_CPG = 20_000  # toy genome-wide CpG count


def simulate_condition(name, mean_cov, true_meth_frac=0.78, conversion=0.9965):
    """Simulate per-CpG calls for one library.

    mean_cov:       target mean coverage (10 ng high, 0.1 ng low + noisy)
    true_meth_frac: fraction of CpGs truly methylated (human ~0.7-0.8)
    conversion:     EM-seq conversion efficiency
    """
    cov = RNG.poisson(mean_cov, N_CPG)
    is_meth = RNG.random(N_CPG) < true_meth_frac
    # methylated C is protected -> read methylated ~= conversion;
    # unmethylated C is converted -> residual failure = 1 - conversion.
    p_meth = np.where(is_meth, conversion, 1 - conversion)
    meth_reads = RNG.binomial(cov, p_meth)
    df = pd.DataFrame({
        "chrom": "chr_toy",
        "pos": np.arange(1, N_CPG + 1) * 100,
        "meth_reads": meth_reads,
        "unmeth_reads": cov - meth_reads,
    })
    df = df[df.meth_reads + df.unmeth_reads > 0].reset_index(drop=True)
    df.to_csv(OUT / f"cpg_{name}.tsv", sep="\t", index=False)
    return df


def simulate_spikein(name, n_sites, meth_target, mean_cov=200):
    cov = RNG.poisson(mean_cov, n_sites)
    meth_reads = RNG.binomial(cov, meth_target)
    df = pd.DataFrame({
        "control": name,
        "pos": np.arange(1, n_sites + 1),
        "meth_reads": meth_reads,
        "unmeth_reads": cov - meth_reads,
    })
    df.to_csv(OUT / f"control_{name}.tsv", sep="\t", index=False)
    return df


if __name__ == "__main__":
    simulate_condition("10ng", mean_cov=30, conversion=0.9968)
    simulate_condition("0p1ng", mean_cov=6, conversion=0.9959)   # ultra-low input
    # unmethylated lambda -> reads ~0% methylated -> conversion efficiency
    simulate_spikein("lambda", n_sites=1500, meth_target=0.003)
    # CpG-methylated pUC19 -> reads ~high% methylated -> protection / over-conversion
    simulate_spikein("puc19", n_sites=1500, meth_target=0.972)
    print("wrote synthetic data to", OUT)
