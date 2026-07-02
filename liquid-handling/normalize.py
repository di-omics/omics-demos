#!/usr/bin/env python3
"""Normalize a 96-well DNA plate to a target mass, then row-pool, on a Hamilton
STARlet driven by PyLabRobot's chatterbox (simulation) backend -- no hardware.

Steps:
  1. Compute per-well sample + diluent volumes to hit TARGET_NG at FINAL_UL, with
     edge handling: wells too dilute to reach the target, and wells so concentrated
     they'd need a sub-microliter transfer, are both flagged.
  2. Execute the transfers column-by-column on an 8-channel head. Diluent (water)
     reuses one tip set; sample transfers get fresh tips per column (no carryover).
  3. Row-pool a fixed volume from every normalized well into 8 pooled fractions.

Artifacts: transfer worklist (CSV), a per-well report (TSV), and the full deck
action log. Swap LiquidHandlerChatterboxBackend for the STAR backend to run for real.
"""
import asyncio, contextlib, io
import numpy as np, pandas as pd
from pathlib import Path

from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends.chatterbox import LiquidHandlerChatterboxBackend
from pylabrobot.resources import (
    STARLetDeck, TIP_CAR_288_A00, PLT_CAR_L5MD_A00,
    hamilton_96_tiprack_300uL_filter, Cor_96_wellplate_360ul_Fb,
)

DATA = Path(__file__).parent / "data"
TARGET_NG, FINAL_UL = 50.0, 25.0
MIN_TX, MIN_DILUENT = 1.0, 1.0
MAX_SAMPLE = FINAL_UL - MIN_DILUENT
POOL_UL = 5.0
ROWS, COLS = list("ABCDEFGH"), list(range(1, 13))


def build_plan(df):
    conc = df.conc_ng_ul.values
    needed = TARGET_NG / conc
    sample = np.clip(needed, MIN_TX, MAX_SAMPLE)
    status = np.where(needed > MAX_SAMPLE, "dilute",
             np.where(needed < MIN_TX, "min_vol", "ok"))
    diluent = np.round(FINAL_UL - sample, 1)
    sample = np.round(sample, 1)
    df = df.copy()
    df["sample_ul"], df["diluent_ul"] = sample, diluent
    df["out_mass_ng"] = (conc * sample).round(1)
    df["status"] = status
    return df


def worklist(df):
    rows = []
    for _, r in df.iterrows():
        rows.append(("normalize_diluent", "diluent", r.well, "normalized", r.well, r.diluent_ul))
        rows.append(("normalize_sample", "samples", r.well, "normalized", r.well, r.sample_ul))
        rows.append(("row_pool", "normalized", r.well, "pool", f"{r.row}1", POOL_UL))
    return pd.DataFrame(rows, columns=[
        "step", "source_labware", "source_well", "dest_labware", "dest_well", "volume_ul"])


async def run_deck(df):
    deck = STARLetDeck()
    lh = LiquidHandler(backend=LiquidHandlerChatterboxBackend(num_channels=8), deck=deck)

    tip_car = TIP_CAR_288_A00(name="tip_carrier")
    tip_car[0] = tips_dil = hamilton_96_tiprack_300uL_filter(name="tips_diluent")
    tip_car[1] = tips_smp = hamilton_96_tiprack_300uL_filter(name="tips_sample")
    tip_car[2] = tips_pool = hamilton_96_tiprack_300uL_filter(name="tips_pool")
    deck.assign_child_resource(tip_car, rails=1)

    plt_car = PLT_CAR_L5MD_A00(name="plate_carrier")
    plt_car[0] = src = Cor_96_wellplate_360ul_Fb(name="samples")
    plt_car[1] = dil = Cor_96_wellplate_360ul_Fb(name="diluent")
    plt_car[2] = dst = Cor_96_wellplate_360ul_Fb(name="normalized")
    plt_car[3] = pool = Cor_96_wellplate_360ul_Fb(name="pool")
    deck.assign_child_resource(plt_car, rails=10)

    by_col = {c: df[df.col == c].sort_values("row") for c in COLS}

    await lh.setup()

    # 1. diluent: reuse one tip set across all columns (water, no carryover)
    await lh.pick_up_tips(tips_dil["A1:H1"])
    for c in COLS:
        v = list(by_col[c].diluent_ul.values)
        await lh.aspirate(dil[f"A{c}:H{c}"], vols=v)
        await lh.dispense(dst[f"A{c}:H{c}"], vols=v)
    await lh.drop_tips(tips_dil["A1:H1"])

    # 2. sample: fresh tips per column
    for c in COLS:
        v = list(by_col[c].sample_ul.values)
        await lh.pick_up_tips(tips_smp[f"A{c}:H{c}"])
        await lh.aspirate(src[f"A{c}:H{c}"], vols=v)
        await lh.dispense(dst[f"A{c}:H{c}"], vols=v)
        await lh.drop_tips(tips_smp[f"A{c}:H{c}"])

    # 3. row-pool: fixed volume from every normalized well into 8 pooled fractions
    for c in COLS:
        v = [POOL_UL] * 8
        await lh.pick_up_tips(tips_pool[f"A{c}:H{c}"])
        await lh.aspirate(dst[f"A{c}:H{c}"], vols=v)
        await lh.dispense(pool["A1:H1"], vols=v)   # accumulate into column 1 (row pools)
        await lh.drop_tips(tips_pool[f"A{c}:H{c}"])


def main():
    df = pd.read_csv(DATA / "concentrations.tsv", sep="\t")
    df = build_plan(df)

    # execute, capturing the (verbose) chatterbox log to a file
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        asyncio.run(run_deck(df))
    log = buf.getvalue()
    (DATA / "deck_actions.log").write_text(log)

    wl = worklist(df)
    df.to_csv(DATA / "normalization_report.tsv", sep="\t", index=False)
    wl.to_csv(DATA / "worklist.csv", index=False)

    # accounting
    n = len(df)
    tips_used = 8 + n + n            # diluent set + one sample tip/well + one pool tip/well
    vol_moved = (df.diluent_ul.sum() + df.sample_ul.sum() + POOL_UL * n)
    counts = df.status.value_counts().to_dict()
    ok = df[df.status == "ok"]

    print("=== Hamilton STARlet normalization + row-pool (synthetic, chatterbox sim) ===")
    print(f"plate: {n} wells, 12 columns, 8-channel head")
    print(f"target {TARGET_NG:.0f} ng in {FINAL_UL:.0f} uL   |   pool {POOL_UL:.0f} uL/well -> 8 row fractions")
    print(f"well status: " + ", ".join(f"{k}={v}" for k, v in sorted(counts.items())))
    if len(ok):
        print(f"on-target wells: output mass {ok.out_mass_ng.mean():.1f} ng, CV {100*ok.out_mass_ng.std()/ok.out_mass_ng.mean():.1f}%")
    print(f"tips used: {tips_used}   |   total volume moved: {vol_moved:.0f} uL   |   deck actions logged: {len(log.splitlines())} lines")

    flagged = df[df.status != "ok"]
    if len(flagged):
        print(f"\nflagged wells ({len(flagged)}):")
        print(flagged[["well", "conc_ng_ul", "sample_ul", "out_mass_ng", "status"]].to_string(index=False))

    print(f"\nwrote {DATA/'normalization_report.tsv'}, {DATA/'worklist.csv'}, {DATA/'deck_actions.log'}")


if __name__ == "__main__":
    main()
