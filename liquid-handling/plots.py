#!/usr/bin/env python3
"""Plot for the normalization demo: input concentration vs normalized output mass."""
import pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

ROOT = Path(__file__).parent
DATA, ASSETS = ROOT / "data", ROOT / "assets"; ASSETS.mkdir(exist_ok=True)
TARGET = 50.0


def main():
    df = pd.read_csv(DATA / "transfer_plan.tsv", sep="\t")
    fig, ax = plt.subplots(1, 2, figsize=(9.5, 3.6))
    ax[0].bar(df.well, df.conc_ng_ul, color="#8a8a8a")
    ax[0].set(title="Input concentration", ylabel="ng/µL", xlabel="well")
    ax[1].bar(df.well, df.out_mass_ng, color="#3b6ea5")
    ax[1].axhline(TARGET, ls="--", color="#a53b6e", lw=1, label=f"target {TARGET:.0f} ng")
    ax[1].set(title="Normalized output mass", ylabel="ng", xlabel="well"); ax[1].legend(frameon=False)
    fig.tight_layout(); out = ASSETS / "normalization_qc.png"; fig.savefig(out, dpi=130)
    print("wrote", out)


if __name__ == "__main__":
    main()
