# Apple Silicon compatibility evidence — 2026-09-21

Scope: dependency definitions under `workflow/workflow/envs/`, checked before installation. This is package/platform evidence; successful solving, executable checks and actual workflow execution must be recorded separately.

| Component | Current availability / relevant constraint | Authoritative evidence |
| --- | --- | --- |
| Snakemake / snakemake-minimal | Platform-independent Conda package; current recipe requires Python >=3.11 | [Bioconda recipe](https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/snakemake/meta.yaml) |
| minimap2 | Native osx-arm64; package page lists 2.31 | [Bioconda package](https://anaconda.org/bioconda/minimap2) |
| samtools | Native osx-arm64; recipe lists 1.24 | [Bioconda recipe](https://bioconda.github.io/recipes/samtools/README.html) |
| pysam | Native osx-arm64; recipe lists 0.24.1; Python 3.12 builds exist | [Bioconda recipe](https://bioconda.github.io/recipes/pysam/README.html), [native package index](https://conda.anaconda.org/bioconda/osx-arm64/) |
| NanoPlot | Platform-independent package, currently 1.48.0; native compiled dependencies must resolve | [Bioconda package](https://anaconda.org/bioconda/nanoplot), [recipe dependencies](https://bioconda.github.io/recipes/nanoplot/README.html) |
| filtlong | Native osx-arm64; 0.3.1 | [Bioconda package](https://anaconda.org/bioconda/filtlong) |
| whatshap | Native osx-arm64; package page lists 2.8 | [Bioconda package](https://anaconda.org/bioconda/whatshap) |
| bcftools / htslib | Native osx-arm64; recipes list 1.24 | [bcftools recipe](https://bioconda.github.io/recipes/bcftools/README.html), [htslib recipe](https://bioconda.github.io/recipes/htslib/README.html) |
| Clair3 | Native osx-arm64 build `clair3-2.0.3-py311h9aa1f4a_0.conda`, uploaded 2026-09-09, approximately 402.5 MB | [native package index](https://conda.anaconda.org/bioconda/osx-arm64/) |

Clair3 2.0.3's recipe explicitly builds for osx-arm64, includes an ARM64 PyPy distribution, and requires Python 3.11. Keep its environment separate from the repository's Python 3.12 BAM QC/report environments. The recipe tests `longphase --version` and `run_clair3.sh -v`; these are not an end-to-end calling test. [Clair3 recipe](https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/clair3/meta.yaml)

Clair3 v2 moved from TensorFlow to PyTorch. Older TensorFlow `.index`/`.data` checkpoints do not work with v2. Compatible checkpoints are `pileup.pt` and `full_alignment.pt`; the package bundles converted models, including the configured `r1041_e82_400bps_hac_v520`, under `${CONDA_PREFIX}/bin/models/`. Chemistry/basecaller matching still requires the sequencing run metadata. FASTQ-remapped BAMs do not preserve Dorado move-table tags, so do not select a `_with_mv` model for this input route. Upstream directs GPU users to its manual setup; native CPU execution is the appropriate first validation path. [Clair3 upstream documentation](https://github.com/HKU-BAL/Clair3#pre-trained-models), [model packaging recipe](https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/clair3/build.sh)

NanoPlot's current recipe requires Plotly >=6.1.1 and Kaleido >=1. Kaleido 1 uses a separate Chrome/Chromium installation for static plots. If a compatible browser cannot run, NanoPlot's documented `--no_static` option retains the HTML/statistics route. Check actual report output before marking QC complete. [Plotly browser requirement](https://plotly.com/python/static-image-export/#chrome), [NanoPlot options](https://github.com/wdecoster/NanoPlot#usage)

Recommended setup sequence: solve explicitly for `osx-arm64` with conda-forge before bioconda and strict channel priority; create workspace-local prefixes without modifying base; record explicit package exports; validate executable architectures/versions and imports; perform workflow dry runs before processing KEL. Native package availability does not establish model suitability or biological validation.
