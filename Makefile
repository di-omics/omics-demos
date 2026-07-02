PY ?= python

.PHONY: all emseq umi species flex liquid chromatin rna clean

all: emseq umi species flex liquid chromatin rna

emseq:
	cd emseq-methylation && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

umi:
	cd umi-dedup && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

species:
	cd species-transcriptome && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

flex:
	cd flex-rna && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

liquid:
	cd liquid-handling && $(PY) generate_data.py && $(PY) normalize.py && $(PY) plots.py

chromatin:
	cd chromatin-browser && $(PY) make_preview.py
	@echo ">> open chromatin-browser/index.html in a browser for the interactive version"

rna:
	cd rna-seq && $(PY) generate_data.py && $(PY) analyze.py && $(PY) plots.py

clean:
	rm -rf */data */__pycache__
