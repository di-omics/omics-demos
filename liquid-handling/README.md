# liquid-handling

Targeted PCR library prep on a Hamilton STAR, driven by **PyLabRobot** in simulation
(chatterbox) mode - the whole protocol executes and logs every deck action with **no
hardware**. This is the automation layer under the rest of the omics stack: hand it a
plate of templates, it runs library prep and hands back an indexed library plus a
worklist.

The volumes and deck layout here are **generic and illustrative**, not a validated
method. All data is synthetic - no real samples.

## Run

```bash
pip install -r requirements.txt   # includes pylabrobot
make liquid
```

## What it does

One column, 8-channel, three phases:

1. **PCR1 master mix** - distribute PCR1 master mix onto the templates.
2. **SPRI cleanup** - add beads, move the bead-bound reaction to the magnet position,
   remove the supernatant, two 80% ethanol washes, elute, and collect the eluate.
3. **PCR2 index master mix** - add the indexing master mix to the cleaned product.

Everything is built from a single transfer plan, so the executed deck actions and the
exported worklist can't drift. Fresh tips are taken for every transfer; p50 tips handle
master mix and elution, p300 tips handle beads, supernatant, and ethanol.

## Example output

```
=== Hamilton STAR amplicon-seq library prep (synthetic, chatterbox sim) ===
mode: DRY / simulated   |   column-1, 8-channel   |   3 phases
          phase  transfers  volume_ul
PCR1 master mix          1         20
   SPRI cleanup          9        819
  PCR2 index MM          1         18

total transfers: 11   |   tips used: p50=32, p300=56   |   volume moved: 6856 uL (x8 channels)
```

To run on a real instrument, swap `LiquidHandlerChatterboxBackend` for the STAR
backend behind typed confirmations - the protocol plan is unchanged.

![Deck map and per-step transfers](assets/libprep_qc.png)

## Files

```
generate_data.py   synthesize an 8-sample sheet (column 1) with i7/i5 indices
run_protocol.py    build the transfer plan + execute it on the STAR (chatterbox)
plots.py           Hamilton STAR deck map + per-step transfer chart
```

## Outputs

```
data/worklist.csv            per-transfer worklist (phase, reagent, source, dest, volume, tips)
data/protocol_summary.tsv    transfers + volume per phase
data/deck_actions.log        every simulated deck action
assets/libprep_qc.png        deck map + per-step transfers
```
