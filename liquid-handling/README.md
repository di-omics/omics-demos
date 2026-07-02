# liquid-handling

DNA normalization on a Hamilton STARlet, driven by **PyLabRobot** in simulation
(chatterbox) mode - the whole protocol runs and logs every deck action with **no hardware**.

All data is synthetic - no real samples.

## Run

```bash
# from the repo root
pip install -r requirements.txt   # includes pylabrobot
make liquid
```

## What it does

- Reads per-well DNA concentrations (8 wells).
- Computes sample + diluent volumes to hit a fixed **target mass at a fixed final volume**.
- Lays out a STARlet deck (tip carrier + 3 Corning plates) and executes an 8-channel
  diluent-then-sample transfer. The chatterbox backend prints tips, aspirates, and dispenses.

## Example output

```
well  conc_ng_ul  sample_ul  diluent_ul  out_mass_ng
  A1        49.1        1.0        24.0         49.1
  E1         7.0        7.1        17.9         49.7
  H1         6.5        7.7        17.3         50.1

target 50.0 ng in 25.0 uL | output mass CV: 1.6%
```

Inputs spanning 6.5-49.2 ng/µL collapse to ~50 ng out (CV 1.6%). To run on a real
instrument, swap `LiquidHandlerChatterboxBackend` for the STAR backend - the protocol is unchanged.

![Normalization](assets/normalization_qc.png)

## Files

```
generate_data.py   synthesize per-well concentrations
normalize.py       compute transfer plan + run the STARlet protocol (chatterbox)
plots.py           input concentration vs normalized output mass
```
