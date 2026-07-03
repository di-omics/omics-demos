# demux-index-hopping

Quantify index hopping on a patterned flow cell, in one runnable example.
Simulates dual-indexed reads for 8 samples, injects a controlled hopping rate,
then recovers the hopping percentage from the i7 x i5 count matrix.

All data is synthetic -- no real sequencing data.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make demux
```

## What it does

- Generates ~50,000 dual-indexed reads across 8 samples, each with a unique (i7, i5) pair.
- Injects ~2% index hopping (one correct index, one swapped) and ~1% undetermined reads.
- Assigns reads by exact index-pair match and builds an i7 x i5 count matrix.
- Quantifies the hopping rate (off-diagonal / total classified) and flags unexpected combos.

## Example output

```
Total reads:           50000
Classified (known i7xi5): 49490
  On-diagonal (correct):  48498
  Off-diagonal (hopped):  992
Undetermined:          1502
Estimated hopping rate:   2.00%
```

Hopping rate recovered close to the injected value (~2%).

![Demux QC](assets/demux_qc.png)

## How it works

Build the i7 x i5 count matrix and compute the hopping rate from off-diagonal entries (from `analyze.py`):

```python
matrix = np.zeros((n, n), dtype=int)
for _, r in reads.iterrows():
    ri = i7_idx.get(r.i7)
    ci = i5_idx.get(r.i5)
    if ri is not None and ci is not None:
        matrix[ri, ci] += 1

diag_total = int(np.trace(matrix))
off_diag_total = int(matrix.sum() - diag_total)
hop_rate = off_diag_total / classified if classified else 0
```

## Files

```
generate_data.py   simulate dual-indexed reads with index hopping + undetermined
analyze.py         assign reads, build hopping matrix, quantify hopping rate
plots.py           i7 x i5 heatmap + per-sample bar chart
```
