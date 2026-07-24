#!/usr/bin/env python3
"""PCR enrichment on a Hamilton STAR with the PyLabRobot simulation
backend -- no hardware. Column-1, 8-channel, three phases:

  1. PCR enrichment   -> add master mix to templates
  2. SPRI cleanup     -> beads, move to magnet, remove supernatant, 2x ethanol
                         wash, elute, collect
  3. Indexing PCR     -> add indexing master mix to the cleaned product

The whole run is built from one transfer plan, so the executed deck actions and
the exported worklist can never drift. Volumes and layout are generic and
illustrative, not a validated method. All data is synthetic. To run for real,
swap the chatterbox backend for the STAR backend behind typed confirmations.
"""
import asyncio, contextlib, io
import pandas as pd
from pathlib import Path

from pylabrobot.liquid_handling import LiquidHandler
from pylabrobot.liquid_handling.backends.chatterbox import LiquidHandlerChatterboxBackend
from pylabrobot.resources import (
    STARDeck, TIP_CAR_288_A00, PLT_CAR_L5MD_A00,
    hamilton_96_tiprack_50uL, hamilton_96_tiprack_300uL_filter, Cor_96_wellplate_360ul_Fb,
)

DATA = Path(__file__).parent / "data"
DRY_RUN = True   # chatterbox = simulated/dry; real runs use the STAR backend + confirmations

# generic, illustrative volumes (uL)
ENRICHMENT_MM, INDEX_MM, SPRI, BEAD_MOVE, SUP, ETOH, ELUTE = 20, 18, 45, 68, 66, 150, 20

# column-1 ranges by labware column
C = {1: "A1:H1", 2: "A2:H2", 3: "A3:H3", 4: "A4:H4"}

# transfer plan: single source of truth for execution + worklist.
# (phase, action, reagent, src_labware, src_range, dst_labware, dst_range, volume, tips)
PLAN = [
    ("PCR enrichment", "transfer", "PCR enrichment mix", "reagent", C[1], "work", C[1], ENRICHMENT_MM, "p50"),
    ("PCR enrichment", "note", "thermocycle PCR enrichment (offline)", None, None, None, None, None, None),

    ("SPRI cleanup", "transfer", "SPRI beads 1.8x", "reservoir", C[1], "work", C[1], SPRI, "p300"),
    ("SPRI cleanup", "note", "mix, bind 5 min", None, None, None, None, None, None),
    ("SPRI cleanup", "transfer", "bead-bound reaction", "work", C[1], "magnet", C[1], BEAD_MOVE, "p300"),
    ("SPRI cleanup", "note", "engage magnet, 3 min", None, None, None, None, None, None),
    ("SPRI cleanup", "transfer", "supernatant", "magnet", C[1], "waste", C[1], SUP, "p300"),
    ("SPRI cleanup", "transfer", "80% ethanol wash 1", "reservoir", C[2], "magnet", C[1], ETOH, "p300"),
    ("SPRI cleanup", "transfer", "ethanol 1", "magnet", C[1], "waste", C[1], ETOH, "p300"),
    ("SPRI cleanup", "transfer", "80% ethanol wash 2", "reservoir", C[3], "magnet", C[1], ETOH, "p300"),
    ("SPRI cleanup", "transfer", "ethanol 2", "magnet", C[1], "waste", C[1], ETOH, "p300"),
    ("SPRI cleanup", "note", "air dry 2 min, disengage magnet", None, None, None, None, None, None),
    ("SPRI cleanup", "transfer", "elution buffer", "reservoir", C[4], "magnet", C[1], ELUTE, "p50"),
    ("SPRI cleanup", "note", "mix, elute, re-engage magnet", None, None, None, None, None, None),
    ("SPRI cleanup", "transfer", "eluate", "magnet", C[1], "work", C[2], ELUTE, "p50"),

    ("Indexing PCR", "transfer", "Indexing PCR mix", "reagent", C[2], "work", C[2], INDEX_MM, "p50"),
    ("Indexing PCR", "note", "thermocycle indexing PCR (offline)", None, None, None, None, None, None),
]


def worklist():
    rows = [(p, r, sl, sr, dl, dr, v, t) for (p, a, r, sl, sr, dl, dr, v, t) in PLAN if a == "transfer"]
    return pd.DataFrame(rows, columns=[
        "phase", "reagent", "source", "source_wells", "dest", "dest_wells", "volume_ul", "tips"])


async def run():
    deck = STARDeck()
    lh = LiquidHandler(backend=LiquidHandlerChatterboxBackend(num_channels=8), deck=deck)

    tip_car = TIP_CAR_288_A00(name="tip_carrier")
    tip_car[0] = rack_p50 = hamilton_96_tiprack_50uL(name="p50_tips")
    tip_car[1] = rack_p300 = hamilton_96_tiprack_300uL_filter(name="p300_tips")
    deck.assign_child_resource(tip_car, rails=1)

    plt_car = PLT_CAR_L5MD_A00(name="plate_carrier")
    plt_car[0] = reagent = Cor_96_wellplate_360ul_Fb(name="reagent")
    plt_car[1] = work = Cor_96_wellplate_360ul_Fb(name="work")
    plt_car[2] = magnet = Cor_96_wellplate_360ul_Fb(name="magnet")
    plt_car[3] = reservoir = Cor_96_wellplate_360ul_Fb(name="reservoir")
    plt_car[4] = waste = Cor_96_wellplate_360ul_Fb(name="waste")
    deck.assign_child_resource(plt_car, rails=10)

    labware = {"reagent": reagent, "work": work, "magnet": magnet,
               "reservoir": reservoir, "waste": waste}
    racks = {"p50": rack_p50, "p300": rack_p300}
    tip_col = {"p50": 0, "p300": 0}

    await lh.setup()
    phase = None
    for (ph, action, label, sl, sr, dl, dr, vol, tip) in PLAN:
        if ph != phase:
            phase = ph; print(f"\n===== {ph} =====")
        if action == "note":
            print(f"# {label}")
            continue
        tip_col[tip] += 1
        col = f"A{tip_col[tip]}:H{tip_col[tip]}"
        print(f"\n-- {label}: {sl} {sr} -> {dl} {dr} ({vol} uL, {tip}) --")
        await lh.pick_up_tips(racks[tip][col])
        await lh.aspirate(labware[sl][sr], vols=[vol] * 8)
        await lh.dispense(labware[dl][dr], vols=[vol] * 8)
        await lh.drop_tips(racks[tip][col])
    return tip_col


def main():
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        tips = asyncio.run(run())
    log = buf.getvalue()
    (DATA / "deck_actions.log").write_text(log)

    wl = worklist()
    wl.to_csv(DATA / "worklist.csv", index=False)
    summ = wl.groupby("phase", sort=False).agg(
        transfers=("volume_ul", "size"), volume_ul=("volume_ul", "sum")).reset_index()
    summ.to_csv(DATA / "protocol_summary.tsv", sep="\t", index=False)

    print("=== Hamilton STAR PCR enrichment (synthetic simulation) ===")
    print(f"mode: {'DRY / simulated' if DRY_RUN else 'REAL'}   |   column-1, 8-channel   |   3 phases")
    print(summ.to_string(index=False))
    print(f"\ntotal transfers: {len(wl)}   |   tips used: p50={tips['p50']*8}, p300={tips['p300']*8}"
          f"   |   volume moved: {wl.volume_ul.sum()*8:.0f} uL (x8 channels)")
    print(f"deck actions logged: {len(log.splitlines())} lines")

    # --- validation: worklist vs plan ---
    expected_transfers = sum(1 for _, a, *_ in PLAN if a == "transfer")
    expected_vol = sum(v for _, a, _, _, _, _, _, v, _ in PLAN if a == "transfer")
    ok_count = len(wl) == expected_transfers
    ok_vol = int(wl.volume_ul.sum()) == expected_vol
    ok_tips = tips["p50"] > 0 and tips["p300"] > 0
    all_ok = ok_count and ok_vol and ok_tips
    print(f"\nValidation vs plan: transfers {len(wl)}/{expected_transfers}, "
          f"volume {int(wl.volume_ul.sum())}/{expected_vol} uL, "
          f"tips OK={ok_tips} -> {'PASS' if all_ok else 'FAIL'}")

    print(f"\nwrote {DATA/'worklist.csv'}, {DATA/'protocol_summary.tsv'}, {DATA/'deck_actions.log'}")


if __name__ == "__main__":
    main()
