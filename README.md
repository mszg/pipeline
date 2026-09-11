# Long-range ONT Snakemake pipeline — v0.4

This workflow is for targeted Oxford Nanopore long-range amplicon sequencing. It currently supports the real ABO test dataset used during development: one patient, one gene, three overlapping PCR amplicons, and many FASTQ chunks per amplicon.

## v0.4 patch

v0.4 makes the primer-defined amplicon coordinates the explicit **variant-calling intervals**. The cleaned BAMs remain alignment-quality filters rather than coordinate-clipped BAMs. This preserves useful alignment context while preventing Clair3 from calling outside the intended PCR targets.

The change follows validation against the real ABO dataset. Primer mapping to `NG_006669.2` gave:

```text
fragment1  NG_006669.2:4228-18285   14,058 bp
fragment2  NG_006669.2:11479-24671  13,193 bp
fragment3  NG_006669.2:19001-32388  13,388 bp
```

These independently defined intervals agreed with the v0.3 high-depth intervals to within 0-3 bp at the amplicon boundaries.

### New in v0.4

- `target_region` is no longer used to discard/retain whole BAM alignments;
- primary/MAPQ-filtered variant BAMs retain their complete alignments;
- long-read phasing BAMs likewise retain alignment context;
- per-amplicon 1-based target coordinates are converted into an auditable target manifest;
- overlapping target regions are merged into a standards-compliant **0-based, half-open BED** for each patient+gene analysis;
- Clair3 receives that BED through `--bed_fn`;
- normalized VCF output is restricted to the same BED again with `bcftools view -R` as a defensive post-calling check;
- `clair3.require_target_regions: true` prevents targeted variant calling from being enabled accidentally when an amplicon lacks a target definition.

v0.3's alignment QC remains unchanged: primary mapping rates, MAPQ distributions, secondary/supplementary counts, length-vs-MAPQ QC, observed core intervals, variant BAMs, and long-read phasing BAMs are all retained.

## Input handling

`fastq_input` may be either a directory or a single FASTQ file. Directories are searched recursively for:

```text
*.fastq
*.fastq.gz
*.fq
*.fq.gz
```

Compressed and uncompressed FASTQs can be mixed. All chunks belonging to one amplicon are streamed into a single staged `.fastq.gz`; the originals are never modified.

## Current workflow

```text
FASTQ chunks
     │
     ├─ combine + NanoPlot QC
     ▼
  minimap2
     │
     ▼
RAW sorted BAM ───────────────────────── retained unchanged
     │
     ├─ alignment QC
     │    ├─ primary mapped %
     │    ├─ MAPQ >=20/30/50
     │    ├─ secondary/supplementary counts
     │    └─ length-vs-MAPQ table
     │
     ├──────────────────────────────┐
     ▼                              ▼
VARIANT BAM                     PHASING BAM
primary + mapped                primary + mapped
MAPQ >=30                       MAPQ >=30
no length cutoff                long-read enriched
full alignment retained         full alignment retained
     │                              │
     └──────────────┐               │
                    ▼               │
primer target regions             │
(samples.tsv)                      │
      │                            │
      ▼                            │
merged target BED                  │
      │                            │
      ▼                            │
   Clair3 --bed_fn                 │
      │                            │
      ▼                            │
 target-restricted VCF ────────────┘
      │                            ▼
      └─────────────────────── WhatsHap
                                   │
                                   ▼
                              phased VCF
```

Clair3/WhatsHap remain disabled by default until the appropriate ONT Clair3 model is configured.

## Current ABO configuration

The primer-derived target coordinates are stored directly in `config/samples.tsv`:

```tsv
sample  gene  amplicon   target_region                  phasing_min_length  core_depth_threshold
P001    ABO   fragment1  NG_006669.2:4228-18285        8000                50
P001    ABO   fragment2  NG_006669.2:11479-24671       8000                50
P001    ABO   fragment3  NG_006669.2:19001-32388       8000                20
```

The full file also contains each FASTQ input and reference path. Coordinates in `samples.tsv` are **1-based inclusive**, matching conventional genomic-region notation. The workflow converts them to BED coordinates automatically.

For the current ABO dataset, the three overlapping target regions merge to:

```text
NG_006669.2  4227  32388
```

in BED format (0-based start, half-open end). The original three amplicon intervals remain recorded in `results/targets/P001__ABO/P001__ABO.target_regions.tsv`.

The per-amplicon `core_depth_threshold` values of 50/50/20 are QC thresholds selected for the current ABO test data; they are not universal biological cutoffs.

## Mapping filters

Defaults are in `config/config.yaml`:

```yaml
mapping:
  preset: "map-ont"
  threads: 8
  variant_min_mapq: 30
  phasing_min_mapq: 30
  phasing_min_length: 8000

qc:
  core_depth_threshold: 20
```

### Variant BAM

The variant BAM excludes:

- unmapped reads;
- secondary alignments;
- supplementary alignments;
- primary alignments below the configured MAPQ threshold.

It deliberately does **not** apply a minimum read-length cutoff by default, because shorter reads can still contain valid local SNP/indel evidence. It also does not clip or restrict alignments to the target coordinates; target restriction is performed explicitly at the variant-calling stage.

### Phasing BAM

The phasing BAM uses the same alignment cleanup and additionally requires a configurable minimum read length. This is intended to enrich for reads carrying long-range linkage information. It does not replace the raw BAM.

## Important output files

```text
results/summary/input_manifest.tsv
results/summary/amplicon_summary.tsv
results/summary/gene_summary.tsv

results/targets/<analysis>/<analysis>.bed
results/targets/<analysis>/<analysis>.target_regions.tsv

results/qc/raw/<unit>/NanoStats.txt
results/qc/alignment/<unit>/<unit>.alignment_qc.tsv
results/qc/alignment/<unit>/<unit>.length_mapq.tsv

results/mapping/amplicons/<unit>/<unit>.sorted.bam
results/mapping/amplicons/<unit>/<unit>.variant.bam
results/mapping/amplicons/<unit>/<unit>.phasing.bam
results/mapping/amplicons/<unit>/<unit>.variant.depth.tsv
results/mapping/amplicons/<unit>/<unit>.core_intervals.tsv

results/mapping/genes/<analysis>/<analysis>.merged.bam
results/mapping/genes/<analysis>/<analysis>.variant.bam
results/mapping/genes/<analysis>/<analysis>.phasing.bam
```

`amplicon_summary.tsv` now contains the key QC metrics directly, including primary mapping percentages, MAPQ distributions, secondary/supplementary counts, cleaned-BAM coverage, and the longest observed high-depth interval.

## Running the patched workflow

From the repository root:

```bash
conda activate bloodgroup
snakemake -n -p
```

Then execute:

```bash
snakemake --cores 8 --software-deployment-method conda -p
```

If v0.2 has already been run, Snakemake will reuse the existing raw FASTQ, NanoPlot, reference, and raw BAM outputs where they are still valid. The new v0.3 rules will create the additional QC and cleaned BAM files.

After the run:

```bash
column -t -s $'\t' results/summary/amplicon_summary.tsv
column -t -s $'\t' results/summary/gene_summary.tsv
```

For the detailed length/MAPQ relationship of fragment 1:

```bash
column -t -s $'\t' \
results/qc/alignment/P001__ABO__fragment1/P001__ABO__fragment1.length_mapq.tsv
```

## Variant calling

Target-aware Clair3 calling is implemented, but disabled by default until `clair3.model_path` points to the appropriate ONT model for the basecalling chemistry/model used for the sequencing run.

When enabled:

```yaml
workflow:
  run_variant_calling: true
  run_phasing: true
  run_consensus: true
```

Clair3 consumes:

```text
results/mapping/genes/P001__ABO/P001__ABO.variant.bam
results/targets/P001__ABO/P001__ABO.bed
```

and is invoked with `--bed_fn`, so variant generation is limited to the primer-defined target union. The normalized VCF is restricted to the same BED once more using `bcftools view -R`.

WhatsHap consumes the resulting target-restricted VCF together with:

```text
results/mapping/genes/P001__ABO/P001__ABO.phasing.bam
```

By default, `clair3.require_target_regions: true`. Therefore variant calling will fail at workflow construction if any configured amplicon lacks `target_region`. Set it to `false` only for an intentionally unrestricted analysis.

## Reference

The current ABO test uses the full RefSeqGene sequence `NG_006669.2` saved as:

```text
resources/references/ABO_reference.fasta
```

The workflow itself is not ABO-specific. Other genes/amplicons can be added as additional rows in `config/samples.tsv`.
