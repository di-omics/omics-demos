# variant-calling

Call SNVs from synthetic short reads and see how sensitivity drops at low allele
fractions. High precision; recall drops at low allele fractions -- the honest outcome.

All data is synthetic -- no real sequencing data.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make variant
```

## What it does

- Generates a ~10,000 bp reference and plants ~30 SNVs with allele fractions from 0.05 to 1.0.
- Tiles ~200x coverage of 150 bp reads, flipping variant bases at the true AF rate.
- Piles up reads, computes per-position allele frequencies, and calls variants above a threshold.
- Scores calls against truth: reports **precision**, **recall**, and per-AF-bin **sensitivity**.

## Example output

```
True variants:   30
Called variants:  38
TP / FP / FN:    28 / 10 / 2
Precision:       0.737
Recall:          0.933

Sensitivity by allele-fraction bin:
   af_bin  total  detected  sensitivity
0.00-0.05      6         4        0.667
0.05-0.10      5         5        1.000
0.10-0.20      2         2        1.000
0.20-0.35      5         5        1.000
0.35-0.60      6         6        1.000
0.60-1.01      6         6        1.000
```

![Variant QC](assets/variant_qc.png)

## How it works

Pile up reads per position and call where alt fraction exceeds a threshold (from `analyze.py`):

```python
def call_variants(counts, ref_bases):
    calls = []
    for pos, ctr in enumerate(counts):
        depth = sum(ctr.values())
        if depth < MIN_DEPTH:
            continue
        ref_base = ref_bases[pos]
        for base, count in ctr.items():
            if base == ref_base:
                continue
            af = count / depth
            if af > AF_THRESHOLD:
                calls.append((pos, ref_base, base, round(af, 4), depth))
    return pd.DataFrame(calls, columns=["pos", "ref", "alt", "obs_af", "depth"])
```

## Files

```
generate_data.py   synthesize reference, plant SNVs, generate reads
analyze.py         pileup, call variants, score precision / recall by AF bin
plots.py           sensitivity curve + true-vs-observed AF scatter
```
