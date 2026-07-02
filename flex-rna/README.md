# flex-rna

Probe-based, multiplexed RNA profiling QC - cell calling and per-sample demultiplexing
for a 10x Flex-style design.

> **Assumption:** `flex-demos` here is read as **10x Genomics Flex**. If it means something
> else, swap the data model - the QC pattern is the same.

All data is synthetic - no real counts.

## Run

```bash
# from the repo root
pip install -r requirements.txt
make flex
```

## What it does

- Simulates probe UMI counts for 4 samples multiplexed by probe barcode, plus empty/ambient barcodes.
- Calls cells with a knee-style UMI threshold and scores the call against ground truth.
- Reports per-sample cell count, median UMIs, and median probes detected.

## Example output

```
UMI threshold for cell call: 584
Called cells: 1192   (true cells kept 1191, empties admitted 1, cells lost 9)

Per-sample (called cells):
           cells  median_umi  median_probes
BC01         298      1974.5          190.0
BC02         296      2017.0          190.0
BC03         299      1987.0          189.0
BC04         299      2146.0          192.0
```

![Flex QC](assets/flex_qc.png)

## Files

```
generate_data.py   simulate multiplexed probe counts + empty barcodes
analyze.py         knee-based cell calling, demux, per-sample QC
plots.py           barcode-rank (knee) + probes-per-cell by sample
```
