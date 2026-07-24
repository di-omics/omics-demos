#!/usr/bin/env python3
"""Generate synthetic methylation-sequencing data for the demo.

Produces per-site methylation calls in three contexts (CpG, CHG, CHH) for two
input conditions (10 ng and 0.1 ng) plus functional unmethylated and methylated
controls. Also generates per-read-position methylation rates for M-bias plotting.

All data is synthetic. No real sequencing data is used.
"""
import numpy as np
import pandas as pd
from pathlib import Path

RNG = np.random.default_rng(0)
OUT = Path(__file__).parent / "data"
OUT.mkdir(exist_ok=True)

N_CPG = 20_000
N_CHG = 15_000
N_CHH = 40_000
READ_LEN = 150


def simulate_condition(name, mean_cov, true_meth_frac=0.78, conversion=0.9965):
    """Simulate per-CpG calls for one library."""
    cov = RNG.poisson(mean_cov, N_CPG)
    is_meth = RNG.random(N_CPG) < true_meth_frac
    p_meth = np.where(is_meth, conversion, 1 - conversion)
    meth_reads = RNG.binomial(cov, p_meth)
    df = pd.DataFrame({
        "chrom": "chr_toy", "pos": np.arange(1, N_CPG + 1) * 100,
        "context": "CpG",
        "meth_reads": meth_reads, "unmeth_reads": cov - meth_reads,
    })
    df = df[df.meth_reads + df.unmeth_reads > 0].reset_index(drop=True)

    # CHG context: very low methylation in mammals (~0.5-2%)
    chg_cov = RNG.poisson(mean_cov, N_CHG)
    chg_meth = RNG.binomial(chg_cov, 1 - conversion)  # near-zero, just conversion failures
    chg = pd.DataFrame({
        "chrom": "chr_toy", "pos": np.arange(1, N_CHG + 1) * 150 + 50,
        "context": "CHG",
        "meth_reads": chg_meth, "unmeth_reads": chg_cov - chg_meth,
    })
    chg = chg[chg.meth_reads + chg.unmeth_reads > 0].reset_index(drop=True)

    # CHH context: also near-zero
    chh_cov = RNG.poisson(mean_cov, N_CHH)
    chh_meth = RNG.binomial(chh_cov, 1 - conversion)
    chh = pd.DataFrame({
        "chrom": "chr_toy", "pos": np.arange(1, N_CHH + 1) * 60 + 25,
        "context": "CHH",
        "meth_reads": chh_meth, "unmeth_reads": chh_cov - chh_meth,
    })
    chh = chh[chh.meth_reads + chh.unmeth_reads > 0].reset_index(drop=True)

    all_ctx = pd.concat([df, chg, chh], ignore_index=True)
    all_ctx.to_csv(OUT / f"meth_{name}.tsv", sep="\t", index=False)

    # also save CpG-only for backward compat
    df.to_csv(OUT / f"cpg_{name}.tsv", sep="\t", index=False)
    return all_ctx


def simulate_spikein(name, n_sites, meth_target, mean_cov=200):
    cov = RNG.poisson(mean_cov, n_sites)
    meth_reads = RNG.binomial(cov, meth_target)
    df = pd.DataFrame({
        "control": name, "pos": np.arange(1, n_sites + 1),
        "meth_reads": meth_reads, "unmeth_reads": cov - meth_reads,
    })
    df.to_csv(OUT / f"control_{name}.tsv", sep="\t", index=False)
    return df


def simulate_mbias(name, read_len=READ_LEN, n_reads=50_000, base_rate=0.78,
                   conversion=0.9965):
    """Simulate per-read-position methylation rates (M-bias).

    Real M-bias shows edge artifacts from end-repair / adapter contamination.
    We model a small uptick at the first and last ~10 bp.
    """
    positions = np.arange(1, read_len + 1)
    # base methylation rate across the read
    rate = np.full(read_len, base_rate * conversion)
    # edge artifact: slight increase at ends (fill-in artifact)
    edge = 10
    rate[:edge] += np.linspace(0.04, 0.0, edge)
    rate[-edge:] += np.linspace(0.0, 0.035, edge)
    rate = np.clip(rate, 0, 1)

    meth_calls = RNG.binomial(n_reads, rate)
    total_calls = np.full(read_len, n_reads)

    df = pd.DataFrame({
        "read_position": positions,
        "meth_calls": meth_calls,
        "total_calls": total_calls,
        "meth_rate": (meth_calls / total_calls).round(4),
    })
    df.to_csv(OUT / f"mbias_{name}.tsv", sep="\t", index=False)
    return df


if __name__ == "__main__":
    simulate_condition("10ng", mean_cov=30, conversion=0.9968)
    simulate_condition("0p1ng", mean_cov=6, conversion=0.9959)
    simulate_spikein("unmethylated", n_sites=1500, meth_target=0.003)
    simulate_spikein("methylated", n_sites=1500, meth_target=0.972)
    simulate_mbias("10ng", conversion=0.9968)
    simulate_mbias("0p1ng", conversion=0.9959, n_reads=10_000)
    print("wrote synthetic data to", OUT)
