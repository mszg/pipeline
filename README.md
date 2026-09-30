# Long-range ONT Snakemake pipeline — v0.7

This Snakemake workflow processes Oxford Nanopore long-range amplicon FASTQ
files into alignment/QC summaries, target-aware small-variant calls, phased
variants, exact haplotype allele support and masked reference-guided haplotype
consensuses. Its final step produces a five-sheet Excel evaluation workbook and
a JSON evidence sidecar. ABO and KEL have been exercised with the configured
RefSeqGene references and targets. Additional genes require their own manifest
rows, reference sequence and intended amplicon coordinates; there is no
implemented allele-name or phenotype interpretation.

**v0.7 adds integrated evaluation reporting.** Variant calling, filtering, phase
guard, haplotagging, allele support and consensus behavior are unchanged by the
reporting feature. Reporting implementation and artifact checks passed against
existing KEL and ABO evidence; this does not establish biological genotype/phase
truth or validated assay acceptance thresholds. Current KEL/ABO QC is
**NOT ASSESSED** because the required minimum target-breadth fraction is
unconfigured. See [release notes](docs/releases/v0.7.md) for the complete version
scope, checks rerun for publication and limitations.

The source distribution contains code, environment specifications and synthetic
fixtures. Supply local sequencing inputs/references and a private sample
manifest; checked-in sample tables are configuration examples. Real sequencing
files, generated outputs/reports, local environments, logs and sample-result
evidence are excluded from this release tree. Previously tracked evidence was
retained locally; existing Git history and tags were not rewritten.

## Current workflow

```text
FASTQ chunks
     |
     +-- combine + NanoPlot QC
     v
  minimap2
     |
     v
RAW sorted BAM ------------------------- retained unchanged
     |
     +-- alignment QC
     |
     +---------------------------+
     v                           v
VARIANT BAM                  PHASING BAM
primary + mapped            primary + mapped
MAPQ >= threshold           MAPQ >= threshold
no length cutoff            long-read enriched
     |                           |
     v                           |
primer-defined target BED       |
     |                           |
     v                           |
Clair3 --bed_fn                  |
     |                           |
     v                           |
complete normalized VCF         |
(multiallelic preserved)        |
     |                           |
     +-- reporting/catalogue    |
     |                           |
     v                           |
biallelic phasing-ready VCF ----+
     |
     v
  WhatsHap phase
     |
     +-- phased VCF
     +-- phasing QC
     |
     v
  WhatsHap haplotag
     |
     +-- HP1 BAM / HP2 BAM
     +-- haplotag QC
     +-- haplotype coverage QC
     |
     v
haplotype-specific exact allele support
     |
     +-- variant_support.tsv
     +-- uncertain_variants.tsv
     +-- consensus-ready phased VCF
     |
     v
reference-guided HP1 / HP2 consensus
(unsequenced, low-depth and unresolved regions masked with N)
     |
     v
Excel evaluation workbook + JSON evidence
(execution status and configurable technical QC assessed separately)
```

## Automated Excel evaluation report

The default workflow now ends with `evaluation_report`, producing one workbook
per analysis working directory at `results/reports/pipeline_evaluation.xlsx`
and a machine-readable `pipeline_evaluation.json` evidence sidecar. This adds
reporting only: variant calling, filtering, phasing, support decisions and
consensus algorithms and thresholds are unchanged.

| Sheet | Contents |
| --- | --- |
| `Sample_Summary` | Every expected sample/barcode, expected genes, execution status, overall technical QC, each gene's result and all explanations |
| `Gene_QC` | One row per expected sample/gene: input and primary read counts, mapping percentage, intended-target depth/breadth, coverage gaps, amplicon availability, variant/phasing/support counts and target completeness |
| `Variants` | One row per normalized VCF record, preserving multiallelic ALT order, GT, FILTER, QUAL, DP, AD, AF, phased GT/PS and existing support decisions |
| `Haplotypes` | Phase blocks, eligible phased/unphased heterozygotes, HP read support/depth, FASTA length/N counts, reference-coordinate masks, completeness and FASTA links |
| `Run_Info` | Effective configuration, applied rules, definitions, evidence paths/hashes, reference identity, report-code provenance and recorded VCF tool metadata |

Tables have filters, frozen identifiers/headers, wrapped explanations and explicit
colored status text. IDs and external strings are literal text, including leading
zeros and strings beginning with `=`. Missing annotations display `NOT AVAILABLE`;
real zero remains numeric. No long FASTA sequence is embedded. Links point to local
artifacts and need those files to remain available. Excel cell-limit truncation is
explicitly marked, with warnings in `Run_Info`; full evidence remains in JSON and
linked files.

### Generate or regenerate

The strict rule requires all configured analysis targets and every directly
consumed artifact. It is included in `all`. With a configured environment:

```bash
snakemake --cores 2 --use-conda evaluation_report --dry-run
snakemake --cores 2 --use-conda evaluation_report
```

Dependencies are declared in `workflow/envs/evaluation_report.yaml`
(`python=3.12`, `openpyxl=3.1.5`, `pyyaml>=6`). The report itself does not need
Clair3, pysam or BAM rescanning. For Apple Silicon, use the [native setup instructions](deployment/kel-apple-silicon/README.md),
including the reporting-only library installation. From the repository root in
the documented sibling `setup/` and `workflow/` layout:

```bash
# KEL: provide a private manifest pointing to your existing inputs
bash ../setup/run_snakemake.sh evaluation_report \
  --configfile config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml \
  --config samples=/path/to/local-kel.samples.tsv

# ABO, existing results in the isolated ABO working directory
bash ../setup/run_snakemake.sh evaluation_report \
  --directory /path/to/abo-run \
  --configfile /path/to/abo-config/config.yaml \
    /path/to/abo-config/calling.yaml \
    /path/to/abo-config/reconstruction.yaml
```

Add `--dry-run` to inspect the DAG first. To regenerate a current workbook
explicitly, add `--forcerun evaluation_report`. With current prerequisites these
commands execute only reporting; missing/stale prerequisites are handled normally
by Snakemake. `DEPENDENCIES SATISFIED (cached outputs may be reused)` describes
successful dependency resolution, not a claim that old analysis files ran anew.
Malformed or inconsistent required evidence fails the strict report rule.

### Explicit partial report

When an upstream job fails or outputs are missing, the final strict rule cannot
run. This separate command bypasses the DAG, reads the manifest and available
artifacts, and writes `pipeline_evaluation.partial.xlsx` plus JSON:

```bash
../setup/envs/kel-native/bin/python workflow/scripts/evaluation_report.py \
  --partial --workdir . \
  --config config/config.yaml config/kel.yaml config/kel.calling.yaml config/kel.reconstruction.yaml \
    /path/to/local-kel.yaml

../setup/envs/kel-native/bin/python workflow/scripts/evaluation_report.py \
  --partial --workdir /path/to/abo-run \
  --config config/config.yaml /path/to/abo-config/config.yaml \
    /path/to/abo-config/calling.yaml \
    /path/to/abo-config/reconstruction.yaml
```

For the first command, `/path/to/local-kel.yaml` is a private overlay containing
`samples: /path/to/local-kel.samples.tsv`; it must describe the existing run.
Omit that final overlay only when the checked-in example manifest matches your
configured inputs. The standalone command requires `--partial` and **all config layers**, including
`config/config.yaml`; it does not implicitly load the Snakefile. An optional
`--output PATH.xlsx` changes the destination. Reports are marked
`PARTIAL / INCOMPLETE`, retain all manifest-defined samples/genes, identify missing
stages, and expose parsing/integrity errors as `EVIDENCE ERROR`. Available old
files never establish current execution success. A gene can satisfy its configured
criteria on readable evidence, but the partial report's overall status cannot
be PASS. No historical job failure is inferred without a reliable execution
record; missing files alone do not distinguish failed from never-run steps.

### Manifest identity

The existing `sample` column remains the unique technical workflow key;
`sample__gene` and `sample__gene__amplicon` output paths are unchanged. Optional
TSV columns add reporting metadata:

| Column | Meaning |
| --- | --- |
| `sample_id` | Display ID; defaults to `sample`. Repeated `Positive`/`Negative` labels are allowed with distinct technical `sample` IDs. |
| `barcode` | Sample-level barcode as text; leading zeros retained. Never inferred from FASTQ filenames. |
| `run_id` | Run identity as text; falls back to `reporting.run_id`, otherwise unavailable. |
| `sample_type` | `sample` (default), `positive_control` or `negative_control`; not inferred from labels. |

Rows are grouped by run, technical sample, display sample and barcode. Reused
display labels never merge different technical samples. Use distinct technical
sample IDs for distinct samples or reused barcodes in different runs. If one
existing analysis path pools different reporting identities, its gene metrics
are withheld with an ambiguity explanation; changing reporting metadata cannot
unmix existing analysis outputs. Duplicate workflow units remain errors.

The local KEL/ABO evidence used for reporting checks did not contain sample-level
barcodes or sequencing run IDs in its manifests. The reports show this absence. In particular, KEL's two fragment
barcodes describe amplicons of one sample and must not be entered as two separate
sample identities for the already pooled analysis.

### Configure transparent technical QC

No validated gene-wide acceptance cutoffs are available in this repository.
The default required `target_breadth_fraction` criterion therefore has `min: null`
and gives **NOT ASSESSED**. Existing `haplotypes.callable_min_depth: 50` is used
to *measure* target breadth; it does not imply a validated acceptable breadth
fraction or a gene-level pass. Additional measurement depths can be configured
through `reporting.depth_thresholds`; `reporting.assessment_depth` can override
the breadth measurement depth without changing analysis thresholds.

Configure rules in a YAML overlay under
`reporting.qc.sample_types.TYPE.criteria.METRIC`. Gene overrides use
`reporting.qc.genes.GENE.sample_types.TYPE.criteria.METRIC` and override entire
same-named rules. Set a criterion to `null` to remove it. This template installs
no numerical acceptance cutoff:

```yaml
reporting:
  run_id: "0001"                  # quote text identifiers
  depth_thresholds: [20, 50, 100]  # descriptive measurements only
  qc:
    sample_types:
      sample:
        criteria:
          target_breadth_fraction:
            min: null            # supply your justified fraction, in [0, 1]
            required: true
            severity: FAIL
      positive_control:
        criteria: {}             # configure separately for the assay
      negative_control:
        criteria: {}             # does not inherit sample coverage criteria
    genes:
      KEL:
        sample_types:
          sample:
            criteria:
              assigned_read_fraction:
                min: null        # unavailable cutoff => NOT ASSESSED
                required: true
                severity: REVIEW
```

Rules support inclusive `min`/`max`, `equals`, `one_of`, `required`, `severity`
(`FAIL` or `REVIEW`), explanatory `description`/`units`, and
`not_applicable_if: {metric_name: observed_value}`. Metric names are the `Gene_QC`
column names (also in JSON), e.g. `mean_target_depth`, `zero_coverage_bases`,
`variant_records`, `phased_fraction`, `phase_block_count`, `hp1_min_depth`,
`hp2_target_unmasked_fraction` and `breadth_50x_fraction`. Fractions use 0–1;
`mapping_percent` uses 0–100. Rules and observed values are recorded individually.

| Status | Meaning |
| --- | --- |
| PASS | All required applicable configured criteria assessed and passed |
| REVIEW | A configured warning condition requires review, with no higher-precedence condition |
| FAIL | A required configured failure condition failed |
| NOT ASSESSED | Required criteria or evidence unavailable |

Aggregation is **FAIL > NOT ASSESSED > REVIEW > PASS**, across all expected genes.
All reasons survive aggregation. Missing optional evidence is nonblocking;
failure of an optional criterion is REVIEW. `NOT APPLICABLE` criteria are excluded
and an entirely inapplicable/unconfigured policy never becomes PASS. Zero variants
alone is not failure. Without eligible heterozygotes, phase-fraction/block
criteria are not applicable; FASTA existence does not establish phasing. Controls
use only their own configured criteria. In partial control reports, missing
ordinary-sample depth outputs remain visible execution gaps but do not silently
add sample coverage criteria. Parsing/integrity errors still prevent QC PASS.

### Counting and interpretation

- Input counts are raw FASTQ observations from NanoStats. Mapping percentage is
  primary mapped / all primary mapper records, including unmapped, excluding
  secondary/supplementary; FASTQ filtering can change that denominator.
- Retained variant/phasing counts are primary mapped observations in their
  respective gene BAM flagstats. MAPQ applies to both; phasing also applies its
  per-amplicon length cutoff. These BAMs include alignments outside targets.
  Counts are not independent molecules or deduplicated read-ID counts. Repeated
  source paths/inconsistent input provenance cause pooled counts/depth to be
  withheld. Counts are never summed across genes into a sample total.
- Mean/median/minimum depth and breadth use dense per-amplicon `samtools depth
  -aa` outputs over the intended merged target union, including zeros. Overlapping
  target bases count once while coverage from distinct reads is added. Depth is
  aligned-base depth (not deletion/skip depth), minimum BQ 0, with default
  duplicate/QCFAIL exclusions in addition to the primary/MAPQ BAM filters.
  Truncated dense depth is an error, not zero coverage. Zero/low-depth regions
  appear as intervals and counts; exact BED coordinates are in JSON.
- Observed amplicons means expected amplicons with readable alignment QC. It
  does not prove product identity or sufficient yield. Per-amplicon counts and
  settings remain in JSON. Target coordinates must be available from the manifest;
  missing/inconsistent intended targets are not replaced by whole-reference depth.
- VCF/manifest positions are 1-based, with inclusive manifest ends; BED is
  0-based half-open. Each multiallelic record keeps its ordered ALT list, REF+ALT
  AD vector and ALT AF vector. DP need not equal summed AD, AF is not imputed, and
  missing vector elements remain `.`. Eligible phasing counts refer to biallelic
  heterozygotes; additional columns retain nonbiallelic/all-unphased counts.
- Assigned-read fraction uses selected phasing reads as denominator. HP labels
  are local to gene/contig/phase block and do not imply cross-gene or disconnected
  block linkage. No allele dropout is inferred from imbalance.
- Consensus completeness is the fraction of intended **reference positions**
  outside mask BEDs. It is separate from FASTA length/N counts: indels change
  sequence length and outside-target Ns are expected. Supported deletions can
  count as resolved despite low aligned-base depth. ACCEPT/UNRESOLVED remain
  existing support decisions, not truth-set or gene-QC labels.

Report-code commit/hash and environment versions are distinguished from recorded
analysis metadata. The standard outputs lack original analysis commit provenance,
unique-molecule counts and HP median depth; those are not invented. Recorded VCF
headers supply only the tool provenance actually present. There is no implemented
allele-name or phenotype interpretation, independently established genotype/phase
truth, or validated assay acceptance threshold. Technical QC PASS is not biological
or clinical validation.

### Reporting verification

```bash
../setup/envs/kel-native/bin/python -m unittest discover -s test -p 'test_evaluation*.py' -v
```

Fixtures are synthetic and workbooks from those fixtures are labelled
`SYNTHETIC DEMONSTRATION`. Tests cover missing outputs, zero calls/no heterozygotes,
ABO+KEL aggregation, multiallelic annotations, overlap/zero-depth accounting,
reference-coordinate masks, repeated control labels, identifier/formula safety,
control policies, precedence and malformed evidence. Actual KEL/ABO workbook
checks and their limits are recorded in
[Excel reporting verification](docs/validation/Excel_reporting.md).

## Analysis behavior retained in v0.7

The [Apple Silicon setup bundle](deployment/kel-apple-silicon/README.md) contains
source wrappers, native environment locks, input audit/output verifier scripts
and synthetic smoke drivers. Its instructions use the sibling `setup/` and
`workflow/` layout. Locally generated run evidence is not distributed.

Existing WhatsHap haplotagging assigns informative reads to HP1/HP2 using the
phased biallelic VCF. Haplotype coverage and exact REF/ALT support then determine
which candidate alleles can be used in masked reference-guided consensuses.
The existing indel-support corrections and deletion-aware masks are retained.

Reconstruction requires exactly one informative phase set per target contig. Missing phase sets, multiple disconnected blocks, or a target contig without phased heterozygotes fail before haplotagging. Independent contigs have independent HP orientations. The check is recorded in `phase_set_validation.tsv`.

Small normalized insertions/deletions use aligned read sequence across the full allele window, including equivalent repeat placements and an aligned right flank. Wrong anchors, partial reads, ambiguous bases and conflicting sequence count as OTHER. Nearby variation within the window can conservatively make an observation OTHER. Complex alleles and indels without a usable right flank remain unresolved. SNVs use strict pileup base support; placeholders, reference skips and competing indel events are OTHER. The support TSV includes `SUPPORT_METHOD`. The `mpileup_max_depth` setting applies to SNVs; indel counting includes all eligible primary reads in each HP BAM.

The consensus stage does **not** use a hard global QUAL cutoff. During real-data validation, low-QUAL heterozygous SNVs could show reproducible haplotype-specific support, while several nominal homozygous-alt indels showed conflicting repeat-associated evidence. The workflow uses configurable support criteria and classifies variants as `ACCEPT` or `UNRESOLVED`.

Existing configurable support/masking defaults are:

```yaml
haplotypes:
  callable_min_depth: 50
  min_support_depth: 20
  min_het_alt_fraction: 0.30
  min_het_delta: 0.25
  min_hom_alt_fraction: 0.80
  max_other_fraction: 0.25
  require_pass: true
  mpileup_max_depth: 100000
```

These are conservative workflow defaults, not universal biological constants. They are exposed in `config/config.yaml` so that future datasets can be revalidated without changing code.

### Heterozygous variants

For a phased `0|1` or `1|0` call, the ALT allele must be enriched on the haplotype predicted by the phased genotype. The expected ALT haplotype must have an exact ALT fraction of at least `min_het_alt_fraction`, and the difference from the opposite haplotype must be at least `min_het_delta`.

### Homozygous-alt variants

A `1/1` call is accepted only when both HP1 and HP2 independently show strong exact ALT support. This protects the consensus from repeat-associated indels that Clair3 may genotype as homozygous-alt despite substantial competing alignments.

### Multiallelic and unsupported complex sites

Multiallelic records remain in the complete normalized catalogue but are not forced into the consensus. Unsupported complex alleles, non-PASS variants, unphased heterozygotes, and variants with conflicting haplotype-specific support are written to `uncertain_variants.tsv` and masked in the reconstructed sequence.

## Input handling

`fastq_input` may be either a directory or a single FASTQ file. Directories are searched recursively for `*.fastq`, `*.fastq.gz`, `*.fq`, and `*.fq.gz`. All chunks belonging to one amplicon are streamed into a single staged `.fastq.gz`; the originals are never modified.

## Example ABO configuration

The example primer-derived target coordinates are stored in `config/samples.tsv`:

```text
sample  gene  amplicon   target_region                  phasing_min_length  core_depth_threshold
EXAMPLE001    ABO   fragment1  NG_006669.2:4228-18285        8000                50
EXAMPLE001    ABO   fragment2  NG_006669.2:11479-24671       8000                50
EXAMPLE001    ABO   fragment3  NG_006669.2:19001-32388       8000                20
```

The three overlapping products merge to the BED interval `NG_006669.2 4227 32388`.

## Main output files

```text
results/variants/<analysis>/<analysis>.norm.vcf.gz
results/variants/<analysis>/<analysis>.phasing_ready.vcf.gz
results/phasing/<analysis>/<analysis>.phased.vcf.gz

results/haplotypes/<analysis>/<analysis>.haplotagged.bam
results/haplotypes/<analysis>/<analysis>.HP1.bam
results/haplotypes/<analysis>/<analysis>.HP2.bam

results/qc/haplotypes/<analysis>/<analysis>.haplotag_qc.tsv
results/qc/haplotypes/<analysis>/<analysis>.coverage_qc.tsv
results/qc/haplotypes/<analysis>/<analysis>.variant_support.tsv
results/qc/haplotypes/<analysis>/<analysis>.uncertain_variants.tsv
results/qc/haplotypes/<analysis>/<analysis>.uncertain_regions.bed

results/variants/<analysis>/<analysis>.consensus_ready.vcf.gz
results/consensus/<analysis>/<analysis>.haplotype1.fasta
results/consensus/<analysis>/<analysis>.haplotype2.fasta

results/summary/input_manifest.tsv
results/summary/amplicon_summary.tsv
results/summary/gene_summary.tsv
results/reports/pipeline_evaluation.xlsx
results/reports/pipeline_evaluation.json
```

The consensus FASTAs are reference-guided haplotype consensuses, not de novo assemblies. Sequence outside the primer-defined target, positions below the haplotype-specific callable-depth threshold, and unresolved variant spans are masked with `N`. An accepted deletion can exempt its deleted positions from the low-base-depth mask only on its ALT haplotype, with at least `callable_min_depth` exact ALT reads, a callable anchor, and no outside-target or uncertainty mask anywhere in its REF span. This permits supported deletions with zero aligned bases at deleted positions to be applied. For an unresolved insertion, the reference anchor is masked and the complete ambiguity remains documented in `uncertain_variants.tsv`.

Run focused regressions in an environment with Python, pysam, samtools, bcftools and tabix:

```bash
python -m unittest discover -s test -p 'test_*.py' -v
python test/integration_haplotype_consensus.py
```

## Running the workflow

Configure the sample manifest, reference FASTA and matching Clair3 model first.
With Snakemake installed, from the repository root:

```bash
snakemake --cores 2 --software-deployment-method conda --dry-run
snakemake --cores 2 --software-deployment-method conda --printshellcmds
```

For Apple Silicon, use the native setup wrapper described above; a fresh solve
of every rule environment was not tested for this release. Original FASTQs are
read-only inputs. `qc_only` stages input/QC before reference/target setup.

Useful QC views (substitute your technical analysis key):

```bash
column -t -s $'\t' results/qc/phasing/EXAMPLE001__ABO/EXAMPLE001__ABO.phasing_qc.tsv
column -t -s $'\t' results/qc/haplotypes/EXAMPLE001__ABO/EXAMPLE001__ABO.haplotag_qc.tsv
column -t -s $'\t' results/qc/haplotypes/EXAMPLE001__ABO/EXAMPLE001__ABO.coverage_qc.tsv
column -t -s $'\t' results/qc/haplotypes/EXAMPLE001__ABO/EXAMPLE001__ABO.variant_support.tsv | less -S
```

## Reference

The ABO example uses NCBI RefSeqGene `NG_006669.2` stored as `resources/references/ABO_reference.fasta`. Coordinates are RefSeqGene coordinates, not GRCh38 coordinates. The KEL example uses `NG_007492.3` with targets `411-15302` and `14079-28023` (merged `411-28023`). Other genes and amplicons can be added as additional rows in `config/samples.tsv`.
