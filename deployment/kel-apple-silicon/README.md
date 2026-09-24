# KEL Apple Silicon reproducibility bundle

This directory versions setup work originally executed beside the repository.
Deploy it as `setup/` beside a checkout named `workflow/` to use the tested KEL
config paths. It includes code, environment definitions, exact native package
locks and selected validation evidence. Installed environments, package/model
binaries, FASTQs, BAMs, VCFs, consensus FASTAs and full input audits remain local.

`MANIFEST.json` records source and published checksums. Historical machine paths
in evidence use `<WORKSPACE>`. Three formerly absolute script roots now derive
from their deployed location. Numerical results and validation scope are
unchanged. `RUNBOOK.md` contains the original runtime commands; links in deployed
records use the sibling `setup/` and `workflow/` layout.

## Fresh native setup

These instructions require macOS arm64, Git and Conda. The explicit locks target
`osx-arm64` and `noarch`. Read the pre-install assessment in
[apple_silicon_compatibility.md](apple_silicon_compatibility.md).

From a fresh workspace:

```bash
git clone https://github.com/mszg/pipeline.git workflow
python3 -c 'import shutil; shutil.copytree("workflow/deployment/kel-apple-silicon", "setup")'
CONDA_PKGS_DIRS="$PWD/setup/conda-pkgs" /opt/homebrew/bin/conda create -y --prefix "$PWD/setup/envs/kel-native" --file setup/kel-native.osx-arm64.explicit.txt
CONDA_PKGS_DIRS="$PWD/setup/conda-pkgs" /opt/homebrew/bin/conda create -y --prefix "$PWD/setup/envs/clair3-arm64" --file setup/clair3.osx-arm64.explicit.txt
setup/envs/clair3-arm64/bin/python setup/clair3_preflight.py
```

`copytree` refuses an existing `setup` directory. An already validated workspace
needs no redeployment. Core Python 3.12 supplies Snakemake, minimap2, samtools,
bcftools, htslib, pysam, WhatsHap, NanoPlot and filtlong. Separate Clair3 2.0.3 uses
Python 3.11 and bundled PyTorch `pileup.pt`/`full_alignment.pt` HAC v5.2 models.
The preflight checks native executables, model loading, finite CPU inference and
workflow CLI compatibility. The YAML files are readable specifications; the
explicit lock files record the exact packages used, including packaged models.

The wrapper requests two cores, uses the greedy scheduler and incomplete-job
reruns, keeps caches local, and limits BLAS/OpenMP to one thread per process.
It supplies the isolated tools directly. Using Snakemake's
`--software-deployment-method conda` instead requires separate validation of the
rule-specific environments.

## Reference, input audit and staged execution

Retrieve NCBI RefSeqGene `NG_007492.3` into
`workflow/resources/references/KEL_NG_007492.3.fasta`; verify SHA-256
`488d6efe397c2fc5da5ecc65f43453cf6bc2e465d1c2e9beb8a55eaa5fc0074c`.
The exact source URL, all four primer-pair checks and target coordinates are in
[reference_provenance.json](reference_provenance.json). The reference remains a
local resource following the existing repository ignore policy.

The sample table expects `KEL_Barcodes_72_84 3/KEL_Fragment_1` (barcode84) and
`KEL_Barcodes_72_84 3/KEL_Fragment_2` (barcode72) beside `workflow/`. Edit only the
`fastq_input` paths in `workflow/config/kel.samples.tsv` if needed. `KEL` is a
required technical grouping label for one biological sample. Targets are
`411-15302` and `14079-28023`; union 27,613 bp, overlap 1,224 bp. The declared and
observed basecaller is `dna_r10.4.1_e8.2_400bps_hac@v5.2.0`.

From the workspace root:

```bash
setup/envs/kel-native/bin/python setup/audit_kel_inputs.py --source "KEL_Barcodes_72_84 3" --output setup/input_audit
bash setup/run_snakemake.sh qc_only --configfile config/kel.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml --printshellcmds
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --printshellcmds
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --printshellcmds
```

The audit checks full gzip/FASTQ content, hashes, aggregate length/quality/header
metadata and source size/mtime preservation. It emits no read sequences or
individual read identifiers. Keep the full generated audits local. KEL disables
FASTQ filtering/downsampling; MAPQ 30 applies to derived variant/phasing BAMs,
and the additional 8 kb minimum length applies only to phasing.

## Tests and artifact verification

```bash
setup/envs/kel-native/bin/python -m unittest discover -s workflow/test -p 'test_*.py' -v
PATH="$PWD/setup/envs/kel-native/bin:$PATH" setup/envs/kel-native/bin/python workflow/test/integration_haplotype_consensus.py
setup/envs/kel-native/bin/python setup/verify_kel_outputs.py
setup/envs/kel-native/bin/python setup/verify_kel_calling.py
setup/envs/kel-native/bin/python setup/indel_fix/verify_kel_reconstruction.py
```

The real-output verifiers require a completed run and locally generated input
audit; expected KEL counts are dataset-specific. The synthetic integration test
cleans temporary artifacts by default; `--output /path/to/empty/directory` retains
them. Optional core and Clair3 smoke drivers are included:

```bash
setup/envs/kel-native/bin/python setup/smoke/generate_fixture.py
bash setup/run_snakemake.sh all --directory "$PWD/setup/smoke/run" --configfile "$PWD/setup/smoke/config.yaml" --config samples="$PWD/setup/smoke/samples.tsv" --dry-run
bash setup/run_snakemake.sh all --directory "$PWD/setup/smoke/run" --configfile "$PWD/setup/smoke/config.yaml" --config samples="$PWD/setup/smoke/samples.tsv" --printshellcmds
setup/envs/clair3-arm64/bin/python setup/clair3_smoke/run.py
setup/envs/clair3-arm64/bin/python setup/clair3_candidate_smoke/run.py
```

The core fixture has 48 deterministic reads across two overlapping amplicons.
The first Clair3 smoke checks indexed empty-output execution; the second injects
three SNVs and checks exact alleles/genotypes through pileup and full alignment
inference. The scripts write derived fixture/run directories. Historical summary
JSON files accompany them.

Deeper diagnostics can be rerun after KEL reconstruction:

```bash
PATH="$PWD/setup/envs/kel-native/bin:$PATH" setup/envs/kel-native/bin/python setup/indel_fix/compare_legacy_kel.py
setup/envs/kel-native/bin/python setup/indel_fix/leave_one_out.py
```

The legacy comparison runs archived defective functions against identical HP
BAMs in scratch outputs. The second diagnostic withholds 560 CCT>C and 714 G>C
separately from tagging, evaluating support using the original phased genotype.
Both remained accepted. Site 727 was sensitive to tagging changes and retains
its production mask. `indel_fix/before/` is historical test material; active
workflow scripts live under `workflow/workflow/scripts/` in the workspace.

[VALIDATION.md](VALIDATION.md) records 46 successful core jobs, seven calling/
phasing jobs, 11 reconstruction jobs, the final mask rerun, 38 unit tests, exact
synthetic consensus checks and leave-one-out diagnostics. The final dry run
found all outputs current. Of 45 catalogue records, 33 were accepted and 12
masked as unresolved. Both complete haplotype sequences passed independent
reconstruction checks; all 36 original FASTQ hashes remained unchanged.

This establishes computational behavior and internal consistency. Independent
variant/phase truth, switch-error rate, genome-wide mapping specificity and
representative assay sensitivity/specificity remain outside this evidence.
