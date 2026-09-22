# KEL local setup

Repository clone: `../workflow`, upstream `https://github.com/mszg/pipeline`,
initial commit `ec722ae` (v0.6). Raw reads remain in `../KEL_Barcodes_72_84 3`.
One biological sample was confirmed by the user. `KEL` in the required `sample`
column is a technical grouping label; no patient identifier was assigned.

The setup uses native Apple Silicon packages from conda-forge and Bioconda.
`envs/kel-native` contains the core workflow tools and Snakemake;
`envs/clair3-arm64` contains Clair3 2.0.3 with Python 3.11 and PyTorch models.
The existing base/nanopore environments and Homebrew packages were not modified.
The two `*.explicit.txt` files record exact packages for this platform;
`tool_versions.json` and `clair3_preflight.json` record runtime checks.

From the workspace root, run:

```bash
# Reference-independent QC only
bash setup/run_snakemake.sh qc_only --configfile config/kel.yaml --dry-run

# Mapping, filtering, coverage, and summaries
bash setup/run_snakemake.sh all --configfile config/kel.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml --printshellcmds

# Target-aware calling and phasing, using the separate native Clair3 environment
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml --printshellcmds

# Repaired haplotype reconstruction and consensus
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --dry-run
bash setup/run_snakemake.sh all --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml --printshellcmds

# Focused regressions and independent KEL artifact verification
setup/envs/kel-native/bin/python -m unittest discover -s workflow/test -p 'test_*.py' -v
PATH="$PWD/setup/envs/kel-native/bin:$PATH" setup/envs/kel-native/bin/python workflow/test/integration_haplotype_consensus.py
setup/envs/kel-native/bin/python setup/indel_fix/verify_kel_reconstruction.py
```

The wrapper enters the repository, limits requested cores to two, uses workspace
caches, and puts the isolated core tools on PATH. This local setup intentionally
runs without Snakemake's `--software-deployment-method conda`: the original rule
environment files were inspected, but exact local environments supply their
tools. Clair3 runs through `run_clair3.sh` to isolate its Python/runtime from the
core environment. Do not add `--software-deployment-method conda` to these commands
without separately validating the resulting rule environments.

Reference: NCBI RefSeqGene `NG_007492.3`, 28,313 bases, downloaded into
`../workflow/resources/references/KEL_NG_007492.3.fasta`. The full accession/version
is the contig name. These are RefSeqGene coordinates, not GRCh38 coordinates.
All four user-supplied primer pairs exactly match the reference; evidence, source
URL, checksum, and boundaries are in `reference_provenance.json`.

Targets: fragment1 `411-15302`; fragment2 `14079-28023`; union `411-28023`,
27,613 bp; overlap 1,224 bp. The generated BED must be `NG_007492.3 410 28023`.
Basecaller: `dna_r10.4.1_e8.2_400bps_hac@v5.2.0` in every FASTQ header, matched to
the bundled `r1041_e82_400bps_hac_v520` model. Native Clair3 v2 uses converted
PyTorch models; historical ABO numerical results cannot establish equivalence.

`input_audit/` contains full gzip/FASTQ validation, counts, length distributions,
and compressed-original SHA-256 hashes. Fragment2's supplied filenames have gaps
at indices 2 and 5; verify export completeness separately. No FASTQ filtering or
downsampling is enabled. The 8,000 bp cutoff applies only to the derived phasing
BAM, and MAPQ 30 to the variant/phasing BAMs. Raw sorted BAMs are retained.

NanoPlot generates statistics and HTML plots with `--no_static`; browser-driven
static exports are disabled. `smoke/verification.json` concerns a synthetic core
execution test, not biological validation. `verify_kel_outputs.py` checks KEL core
artifacts and rehashes originals against their pre-run hashes.

See `VALIDATION.md` for observed progress and limitations. The core and calling
configs retain staged execution; the additional `kel.reconstruction.yaml` overlay
enables the repaired and revalidated reconstruction. It requires one informative
phase set per target contig, counts indels across complete aligned allele windows,
and correctly handles supported deletions in consensus masks. The 38 regression
tests and known-truth synthetic integration test pass. Archived defective scripts,
before/after comparisons and execution evidence are under `indel_fix/`.

The completed KEL run produced 45 normalized variants. All 41 eligible biallelic
heterozygous calls phased in one block (`PS=516`, positions 516–25592). The final
calling/phasing dry run reports all requested outputs up to date. Reconstruction
assigned 32,100 reads to HP1 and 21,523 to HP2; 32,267 remain unassigned. It accepted
33 variants and masked 12 unresolved candidates. Both complete FASTA sequences
passed independent reconstruction checks; all 36 raw file hashes remain unchanged.
Output FASTAs are in `../workflow/results/consensus/KEL__KEL/`. External variant
and phase truth-set validation remains outstanding.
