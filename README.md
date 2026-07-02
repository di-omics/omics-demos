# omics-demos

Small, self-contained toy examples from an omics + lab-automation stack - sequencing
QC, single-cell/probe RNA, epigenomics, cross-species transcriptomics, and Hamilton
liquid handling. Each demo runs in one command, uses only **synthetic data**, and shows
exactly one idea.

Clean-room by design: no proprietary data or code - every dataset is generated on the fly.

## Quickstart

```bash
pip install -r requirements.txt
make all          # run every demo
# or one at a time:
make emseq        # EM-seq methylation QC
make umi          # UMI deduplication
make species      # cross-species transcriptome
make flex         # 10x Flex cell calling
make liquid       # Hamilton STARlet normalization (PyLabRobot sim)
make chromatin    # chromatin browser preview
make rna          # bulk RNA-seq DE analysis
```

Python 3.10+. `pylabrobot` is only needed for the `liquid` demo.

---

### emseq-methylation - EM-seq methylation QC
Conversion efficiency, CpG protection, coverage, and global methylation across a standard
(10 ng) and ultra-low (0.1 ng) input, using lambda/pUC19 spike-in controls.

![emseq](emseq-methylation/assets/qc_plots.png)

[-> emseq-methylation/](emseq-methylation/)

---

### umi-dedup - UMI deduplication
Recovers true molecule counts from PCR-amplified, error-containing UMI reads with a
1-mismatch directional collapse; reports duplication rate and a saturation curve.

![umi](umi-dedup/assets/umi_qc.png)

[-> umi-dedup/](umi-dedup/)

---

### species-transcriptome - cross-species maternal transcriptome
Human / bovine / marmoset ortholog comparison at zygote and 2-cell stages: cross-species
correlation plus divergent-ortholog ranking.

![species](species-transcriptome/assets/species_qc.png)

[-> species-transcriptome/](species-transcriptome/)

---

### flex-rna - 10x Flex cell calling
Knee-based cell calling and per-sample demultiplexing for a probe-based, multiplexed RNA
design, scored against ground truth.

![flex](flex-rna/assets/flex_qc.png)

[-> flex-rna/](flex-rna/)

---

### liquid-handling - Hamilton STARlet normalization
DNA normalization to a target mass, executed on a STARlet deck via PyLabRobot's simulation
backend - the full protocol runs and logs with no hardware.

![liquid](liquid-handling/assets/normalization_qc.png)

[-> liquid-handling/](liquid-handling/)

---

### chromatin-browser - interactive CUT&Tag track browser
A self-contained genome browser for H3K27me3 and Pol II S5p over a locus, with a live
crosshair readout and called peaks. Opens with no server.

![chromatin](chromatin-browser/assets/preview.png)

[-> chromatin-browser/](chromatin-browser/)

---

### rna-seq - bulk RNA-seq differential expression
CPM normalization, PCA, and per-gene t-test with BH correction on a synthetic
gene x sample count matrix with planted DE genes - measures recovery rate.

![rna-seq](rna-seq/assets/rna_qc.png)

[-> rna-seq/](rna-seq/)

---

## License

MIT - see [LICENSE](LICENSE).
