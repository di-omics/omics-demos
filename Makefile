PY ?= python3

.PHONY: all methylation scrna-umi scrna-cell-calling rna pcr-enrichment chromatin demux cnv variant clean

all: methylation scrna-umi scrna-cell-calling rna pcr-enrichment chromatin demux cnv variant

methylation:
	cd methylation-sequencing && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

scrna-umi:
	cd scrnaseq-umi-dedup && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

scrna-cell-calling:
	cd scrnaseq-cell-calling && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

rna:
	cd rna-seq && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

pcr-enrichment:
	cd pcr-enrichment-automation && $(PY) generate_data.py && $(PY) run_protocol.py && $(PY) plots.py

demux:
	cd demux-index-hopping && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

cnv:
	cd cnv-ploidy && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

variant:
	cd variant-calling && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

chromatin:
	cd chromatin-browser && $(PY) make_preview.py
	@echo ">> open chromatin-browser/index.html in a browser for the interactive version"

clean:
	rm -rf */data */__pycache__
