# Native Apple Silicon setup sources

This directory contains setup scripts, readable environment definitions, native
package locks and synthetic fixtures for macOS arm64. Use it with this release
checked out as `workflow/`, and deploy a copy as sibling `setup/`. The active
pipeline rules and scripts remain in the checkout.

The source release omits historical generated logs, per-sample result documents
and analysis evidence. Local originals and existing Git history are preserved;
their removal from this release tree does not erase earlier publications. Real
FASTQs, BAMs, VCFs, consensus sequences, workbooks, JSON reports, installed
environments and model binaries must remain local. Public reference/primer
metadata is retained in [reference_provenance.json](reference_provenance.json).

## Create the isolated native environments

These commands assume macOS arm64 and Conda at `/opt/homebrew/bin/conda`.
Read the historical platform assessment in
[apple_silicon_compatibility.md](apple_silicon_compatibility.md). From the parent
of the `workflow/` checkout:

```bash
python3 -c 'import shutil; shutil.copytree("workflow/deployment/kel-apple-silicon", "setup")'
CONDA_PKGS_DIRS="$PWD/setup/conda-pkgs" /opt/homebrew/bin/conda create -y --prefix "$PWD/setup/envs/kel-native" --file setup/kel-native.osx-arm64.explicit.txt
CONDA_PKGS_DIRS="$PWD/setup/conda-pkgs" /opt/homebrew/bin/conda create -y --prefix "$PWD/setup/envs/clair3-arm64" --file setup/clair3.osx-arm64.explicit.txt
setup/envs/kel-native/bin/python -m pip install --only-binary=:all: --no-deps 'openpyxl==3.1.5' 'et-xmlfile==2.0.0'
setup/envs/clair3-arm64/bin/python setup/clair3_preflight.py
```

`copytree` refuses an existing `setup/` directory. Do not overwrite an established
workspace; its native environment can receive just the reporting-library install
command above. The explicit locks record the historical core/caller environments
before Excel reporting and deliberately remain unchanged. They do not include
openpyxl or et-xmlfile. The pinned, pure-Python reporting install is an additional
step; PyYAML is already in the core lock. `kel-native.yaml` now also declares
openpyxl for a fresh environment solve, but that updated solve has not been tested.

The core environment contains Python 3.12, Snakemake, minimap2, samtools,
bcftools, pysam, WhatsHap, NanoPlot and filtlong. The separate Clair3 environment
contains Python 3.11, Clair3 2.0.3 and packaged PyTorch HAC v5.2 models. The
preflight checks native executables, model loading, finite CPU inference and the
caller CLI. It writes new local evidence when executed; none is bundled as proof
that a new installation passed.

`run_snakemake.sh` selects the isolated native tools directly, requests two cores,
limits BLAS/OpenMP threads and keeps caches under `setup/`. `run_clair3.sh` selects
the separate caller environment. This wrapper workflow does not automatically
create the rule-specific Conda environments. The dedicated reporting environment
in `workflow/envs/evaluation_report.yaml` is available for normal Snakemake Conda
deployment; a fresh solve/build of that environment was not performed for this
release.

## Configure KEL inputs and targets

Retrieve NCBI RefSeqGene `NG_007492.3` into
`workflow/resources/references/KEL_NG_007492.3.fasta`. The archived FASTA file
SHA-256 is `488d6efe397c2fc5da5ecc65f43453cf6bc2e465d1c2e9beb8a55eaa5fc0074c`.
[reference_provenance.json](reference_provenance.json) contains the public NCBI
retrieval URL, primer sequences and reference-coordinate checks. The reference
is 28,313 bases; primer coordinates are on this accession, not a genome build.

The example `config/kel.samples.tsv` expects input directories
`data/KEL_Fragment_1` and `data/KEL_Fragment_2` inside the checkout. Configure
those paths for the intended dataset without moving or modifying sequencing
originals. Sample metadata and grouping must describe the actual experiment;
example technical labels do not establish biological identity.

The configured conservative target products are `NG_007492.3:411-15302` and
`NG_007492.3:14079-28023`, using 1-based inclusive coordinates. Their union is
27,613 bases, represented by BED `NG_007492.3 410 28023`; the overlap is 1,224
bases. Confirm that these primer-defined products apply to a new experiment.
The example caller model is `r1041_e82_400bps_hac_v520`; match it to the actual
ONT basecalling model.

From the parent workspace, inspect the inputs and dry-run each intended stage
before execution:

```bash
bash setup/run_snakemake.sh qc_only --configfile config/kel.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml --printshellcmds
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --printshellcmds
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --printshellcmds
```

`qc_only` stages inputs and runs raw QC without requiring a reference. `all` now
includes the final Excel/JSON report for the enabled stages. The KEL example
keeps FASTQ filtering disabled, uses MAPQ 30 for variant/phasing BAMs and adds an
8 kb read-length minimum for phasing. The reporting addition does not change
these filters or caller, phase-guard, support or consensus algorithms.

## Reporting existing outputs

With all requested analysis outputs current, this command regenerates only the
report; inspect the dry run first:

```bash
bash setup/run_snakemake.sh evaluation_report --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --dry-run
bash setup/run_snakemake.sh evaluation_report --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --forcerun evaluation_report
```

The report is `workflow/results/reports/pipeline_evaluation.xlsx`, accompanied
by `pipeline_evaluation.json`. When prerequisites are unavailable, explicitly
request an incomplete report instead:

```bash
setup/envs/kel-native/bin/python workflow/workflow/scripts/evaluation_report.py \
  --partial --workdir workflow \
  --config workflow/config/config.yaml workflow/config/kel.yaml \
    workflow/config/kel.calling.yaml workflow/config/kel.reconstruction.yaml
```

Read the main [reporting documentation](../../README.md) before setting QC
criteria. Missing evidence or required thresholds prevents PASS; the existing
per-haplotype depth threshold does not define a validated gene-wide acceptance
fraction. Execution and technical QC are distinct. No biological genotype/phase
truth, named-allele interpretation or phenotype prediction is established here.

## Synthetic checks and diagnostic source

The following tests use synthetic data and temporary outputs. They are separate
from a real-data analysis and do not establish assay sensitivity or specificity:

```bash
setup/envs/kel-native/bin/python -m unittest discover -s workflow/test -p 'test_*.py' -v
PATH="$PWD/setup/envs/kel-native/bin:$PATH" setup/envs/kel-native/bin/python workflow/test/integration_haplotype_consensus.py
setup/envs/kel-native/bin/python setup/smoke/generate_fixture.py
bash setup/run_snakemake.sh all --directory "$PWD/setup/smoke/run" --configfile "$PWD/setup/smoke/config.yaml" --config samples="$PWD/setup/smoke/samples.tsv" --dry-run
bash setup/run_snakemake.sh all --directory "$PWD/setup/smoke/run" --configfile "$PWD/setup/smoke/config.yaml" --config samples="$PWD/setup/smoke/samples.tsv" --printshellcmds
setup/envs/clair3-arm64/bin/python setup/clair3_smoke/run.py
setup/envs/clair3-arm64/bin/python setup/clair3_candidate_smoke/run.py
```

The core smoke fixture generates deterministic reads across two overlapping
amplicons. Caller smoke scripts exercise empty-output handling and injected
variants. The support/consensus integration fixture compares reconstructed
sequences with known synthetic sequences. These commands generate fresh local
outputs; historical JSON summaries and logs are excluded from this source tree.

Additional audit and historical diagnostic scripts remain available as source.
They have dataset-specific assumptions and expected values, so they are not
universal acceptance tests for new samples. `indel_fix/before/` contains archived
counter implementations for regression comparisons; the active algorithms are
under the checkout's `workflow/scripts/`. Keep generated audit metadata and all
real analysis artifacts outside commits and release assets.
