#!/usr/bin/env python3
"""Normalize DNA to a target mass on a Hamilton STARlet, PyLabRobot chatterbox
(simulation) backend -- no hardware required.

Reads per-well concentrations, computes sample + diluent volumes to hit a fixed
target mass at a fixed final volume, then executes the transfers on an 8-channel
head. The chatterbox backend logs every deck action.
"""
import asyncio
import numpy as np, pandas as pd
from pathlib import Path

from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends.chatterbox import LiquidHandlerChatterboxBackend
from pylabrobot.resources import (
    STARLetDeck, TIP_CAR_288_A00, PLT_CAR_L5MD_A00, hamilton_96_tiprack_300uL_filter, Cor_96_wellplate_360ul_Fb,
)

DATA = Path(__file__).parent / "data"
TARGET_NG = 50.0
FINAL_UL = 25.0
MIN_TX, MAX_TX = 1.0, 20.0


def plan(conc):
    sample = np.clip(TARGET_NG / conc, MIN_TX, MAX_TX)
    diluent = np.clip(FINAL_UL - sample, 0, None)
    return sample.round(1), diluent.round(1)


async def main():
    df = pd.read_csv(DATA / "concentrations.tsv", sep="\t")
    sample_v, diluent_v = plan(df.conc_ng_ul.values)
    df["sample_ul"], df["diluent_ul"] = sample_v, diluent_v
    df["out_mass_ng"] = (df.conc_ng_ul * df.sample_ul).round(1)

    deck = STARLetDeck()
    lh = LiquidHandler(backend=LiquidHandlerChatterboxBackend(num_channels=8), deck=deck)

    tip_car = TIP_CAR_288_A00(name="tip_carrier")
    tip_car[0] = tips = hamilton_96_tiprack_300uL_filter(name="tips")
    deck.assign_child_resource(tip_car, rails=1)

    plt_car = PLT_CAR_L5MD_A00(name="plate_carrier")
    plt_car[0] = src = Cor_96_wellplate_360ul_Fb(name="samples")
    plt_car[1] = dil = Cor_96_wellplate_360ul_Fb(name="diluent")
    plt_car[2] = dst = Cor_96_wellplate_360ul_Fb(name="normalized")
    deck.assign_child_resource(plt_car, rails=10)

    await lh.setup()

    col = "A1:H1"
    print("\n--- diluent transfer ---")
    await lh.pick_up_tips(tips[col])
    await lh.aspirate(dil[col], vols=list(diluent_v))
    await lh.dispense(dst[col], vols=list(diluent_v))
    await lh.drop_tips(tips[col])

    print("\n--- sample transfer ---")
    await lh.pick_up_tips(tips["A2:H2"])
    await lh.aspirate(src[col], vols=list(sample_v))
    await lh.dispense(dst[col], vols=list(sample_v))
    await lh.drop_tips(tips["A2:H2"])

    df.to_csv(DATA / "transfer_plan.tsv", sep="\t", index=False)
    print("\n=== normalization plan ===")
    print(df.to_string(index=False))
    print(f"\ntarget {TARGET_NG} ng in {FINAL_UL} uL | output mass CV: "
          f"{100*df.out_mass_ng.std()/df.out_mass_ng.mean():.1f}%")
    print("wrote", DATA / "transfer_plan.tsv")


if __name__ == "__main__":
    asyncio.run(main())
