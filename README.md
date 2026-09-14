# Long-range ONT Snakemake pipeline — v0.5

This workflow is for targeted Oxford Nanopore long-range amplicon sequencing. It is being validated on a real ABO dataset consisting of one patient, one gene, three overlapping long-range PCR amplicons, and multiple FASTQ chunks per amplicon. The workflow itself is not ABO-specific.

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
full alignments retained    full alignments retained
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
  WhatsHap
     |
     +-- phased VCF
     +-- automated phasing QC
```

Haplotype consensus/reconstruction is still disabled by default in v0.5 while handling of unphased, multiallelic, and low-confidence variants is validated.

## v0.5 — validated phasing architecture

v0.5 changes the boundary between variant normalization and phasing.

The complete normalized VCF now **preserves multiallelic records**. A second VCF containing only records with exactly one ALT allele is generated specifically for WhatsHap. This avoids converting a multiallelic `1/2` genotype into duplicate biallelic records at the same coordinate.

During P001 ABO validation, Clair3 produced a multiallelic repeat-associated call at `NG_006669.2:9109`:

```text
REF=CATATATATATATAT
ALT=C,CAT
GT=1/2
QUAL=7.26
```

Splitting this record caused WhatsHap to warn about a duplicate position and skip one representation. Preserving the site in the full catalogue and excluding multiallelic records from the routine phasing input produced a clean biallelic phasing set.

For P001 ABO, the validated result was:

- 67 target-restricted normalized variant records/sites;
- 66 biallelic records suitable for routine WhatsHap phasing;
- 1 preserved multiallelic record;
- 52 usable heterozygous biallelic variants;
- 52/52 heterozygous variants phased (100%);
- one phase set (`PS=9033`);
- one phase block spanning positions 9033-29201;
- 37,820 reads covering the phased variants;
- 37,199 reads covering at least two variants;
- 38 phase-informative reads selected by WhatsHap.

The detailed validation record is stored in `docs/validation/P001_ABO_v0.5.md`.

## v0.4 — target-aware Clair3 calling

v0.4 made the primer-defined amplicon coordinates explicit variant-calling intervals. Cleaned BAMs remain alignment-quality filters rather than coordinate-clipped BAMs, preserving alignment context while preventing Clair3 from calling outside the intended PCR targets.

Primer mapping to `NG_006669.2` for the ABO validation dataset gave:

```text
fragment1  NG_006669.2:4228-18285   14,058 bp
fragment2  NG_006669.2:11479-24671  13,193 bp
fragment3  NG_006669.2:19001-32388  13,388 bp
```

The three target intervals merge to:

```text
NG_006669.2  4227  32388
```

in BED format (0-based start, half-open end). Coordinates in `config/samples.tsv` remain 1-based inclusive.

Clair3 receives the merged BED through `--bed_fn`. The resulting VCF is defensively restricted to the same target with `bcftools view -T` after normalization.

## Input handling

`fastq_input` may be either a directory or a single FASTQ file. Directories are searched recursively for:

```text
*.fastq
*.fastq.gz
*.fq
*.fq.gz
```

Compressed and uncompressed FASTQs can be mixed. All chunks belonging to one amplicon are streamed into a single staged `.fastq.gz`; the originals are never modified.

## Current ABO configuration

The primer-derived target coordinates are stored in `config/samples.tsv`:

```text
sample  gene  amplicon   target_region                  phasing_min_length  core_depth_threshold
P001    ABO   fragment1  NG_006669.2:4228-18285        8000                50
P001    ABO   fragment2  NG_006669.2:11479-24671       8000                50
P001    ABO   fragment3  NG_006669.2:19001-32388       8000                20
```

The full file also contains each FASTQ input and reference path.

## Mapping filters

Defaults are in `config/config.yaml`:

```yaml
mapping:
  preset: "map-ont"
  threads: 8
  variant_min_mapq: 30
  phasing_min_mapq: 30
  phasing_min_length: 8000
```

### Variant BAM

The variant BAM excludes unmapped, secondary, supplementary, and low-MAPQ alignments. It deliberately does not apply a read-length cutoff because shorter reads can still provide useful local SNP/indel evidence. Alignments are not clipped to the target coordinates.

### Phasing BAM

The phasing BAM uses the same alignment cleanup and additionally applies a configurable minimum read length to enrich for reads carrying long-range linkage information.

## Variant calling and phasing

The current validated configuration uses:

```yaml
workflow:
  run_variant_calling: true
  run_phasing: true
  run_consensus: false

clair3:
  platform: "ont"
  model_name: "r1041_e82_400bps_hac_v520"
```

Clair3 consumes the gene-level variant BAM and the primer-derived target BED. v0.5 then produces two normalized variant outputs:

```text
results/variants/<analysis>/<analysis>.norm.vcf.gz
results/variants/<analysis>/<analysis>.phasing_ready.vcf.gz
```

The first is the complete target-restricted normalized catalogue and preserves multiallelic sites. The second contains exactly biallelic records (`bcftools view -m2 -M2`) and is the input to WhatsHap.

WhatsHap uses:

```text
results/variants/<analysis>/<analysis>.phasing_ready.vcf.gz
results/mapping/genes/<analysis>/<analysis>.phasing.bam
results/reference/<analysis>/reference.fasta
```

and produces:

```text
results/phasing/<analysis>/<analysis>.phased.vcf.gz
results/qc/phasing/<analysis>/<analysis>.phasing_qc.tsv
```

The QC table reports variant counts, phased fraction, number of phase sets, and the size and genomic span of the largest phase block.

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
results/qc/phasing/<analysis>/<analysis>.phasing_qc.tsv

results/mapping/amplicons/<unit>/<unit>.sorted.bam
results/mapping/amplicons/<unit>/<unit>.variant.bam
results/mapping/amplicons/<unit>/<unit>.phasing.bam

results/mapping/genes/<analysis>/<analysis>.merged.bam
results/mapping/genes/<analysis>/<analysis>.variant.bam
results/mapping/genes/<analysis>/<analysis>.phasing.bam

results/variants/<analysis>/<analysis>.norm.vcf.gz
results/variants/<analysis>/<analysis>.phasing_ready.vcf.gz
results/phasing/<analysis>/<analysis>.phased.vcf.gz
```

## Running the workflow

From the repository root:

```bash
conda activate bloodgroup
snakemake -n -p
```

Then execute:

```bash
snakemake --cores 8 --software-deployment-method conda -p
```

To inspect the phasing QC result:

```bash
column -t -s $'\t' results/qc/phasing/P001__ABO/P001__ABO.phasing_qc.tsv
```

## Reference

The current ABO validation uses RefSeqGene `NG_006669.2` stored as:

```text
resources/references/ABO_reference.fasta
```

Other genes and amplicons can be added as additional rows in `config/samples.tsv`.
