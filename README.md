# Long-range ONT Snakemake pipeline — v0.6

This workflow is for targeted Oxford Nanopore long-range amplicon sequencing. It is being validated on a real ABO dataset consisting of one patient, one gene, three overlapping long-range PCR amplicons, and multiple FASTQ chunks per amplicon. The workflow itself is not ABO-specific and is intended to remain usable for other targeted long-range sequencing applications.

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
```

## v0.6 — haplotagging and conservative haplotype reconstruction

The [KEL Apple Silicon bundle](deployment/kel-apple-silicon/README.md) contains
tested setup wrappers, native environment locks, audit/verification scripts,
smoke drivers and selected execution evidence. Its instructions retain the
sibling `setup/` and `workflow/` layout used by the KEL configs.
The [complete feature and change description](docs/changes/KEL_apple_silicon_features.md)
lists the implementation changes, retained pipeline capabilities, tests, observed
KEL results and validation limits.

v0.6 adds reconstruction after v0.5 phasing; historical ABO results are summarized below. The KEL code audit found and repaired indel-support and deletion-masking defects; see [KEL revalidation](docs/validation/KEL_indel_revalidation.md) for current executable evidence and limits. WhatsHap `haplotag` assigns phase-informative reads to HP1 or HP2 using the phased biallelic VCF. The workflow then measures target coverage separately for the two haplotypes and evaluates every biallelic candidate using exact REF/ALT support in the HP1 and HP2 read sets.

Reconstruction now requires exactly one informative phase set per target contig. Missing phase sets, multiple disconnected blocks, or a target contig without phased heterozygotes fail before haplotagging. Independent contigs have independent HP orientations. The check is recorded in `phase_set_validation.tsv`.

Small normalized insertions/deletions use aligned read sequence across the full allele window, including equivalent repeat placements and an aligned right flank. Wrong anchors, partial reads, ambiguous bases and conflicting sequence count as OTHER. Nearby variation within the window can conservatively make an observation OTHER. Complex alleles and indels without a usable right flank remain unresolved. SNVs use strict pileup base support; placeholders, reference skips and competing indel events are OTHER. The support TSV includes `SUPPORT_METHOD`. The `mpileup_max_depth` setting applies to SNVs; indel counting includes all eligible primary reads in each HP BAM.

The consensus stage does **not** use a hard global QUAL cutoff. During real-data validation, low-QUAL heterozygous SNVs could show reproducible haplotype-specific support, while several nominal homozygous-alt indels showed conflicting repeat-associated evidence. Instead, v0.6 uses configurable support criteria and classifies variants as `ACCEPT` or `UNRESOLVED`.

Current validation-derived defaults are:

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

Multiallelic records remain in the complete normalized catalogue but are not forced into the consensus in v0.6. Unsupported complex alleles, non-PASS variants, unphased heterozygotes, and variants with conflicting haplotype-specific support are written to `uncertain_variants.tsv` and masked in the reconstructed sequence.

## P001 ABO v0.6 validation

For `P001__ABO`, WhatsHap haplotagging assigned 37,126 of 42,373 phasing reads (87.6%):

```text
HP1          18,412
HP2          18,714
unassigned    5,247
```

Among assigned reads the balance was 49.6% HP1 versus 50.4% HP2. All 37,126 assigned reads carried the validated phase-set ID `PS=9033`.

Across the complete primer-defined target union `NG_006669.2:4228-32388` (28,161 bp), haplotype-specific coverage was complete:

```text
HP1 mean depth   7,715.68x   minimum 71x
HP2 mean depth   7,727.58x   minimum 527x
```

Both haplotypes therefore had 100% target callability at the v0.6 validation threshold of 50x.

Four candidates were retained as unresolved rather than forced into the consensus:

```text
9109   multiallelic 1/2 repeat-associated indel
13629  CA>C         1/1 with conflicting exact allele support
13963  CT>C         1/1, LowQual, conflicting repeat-associated support
25292  G>GACATACAC  1/1 with poor exact ALT support and many competing indels
```

The detailed validation evidence is recorded in `docs/validation/P001_ABO_v0.6.md`.

## Input handling

`fastq_input` may be either a directory or a single FASTQ file. Directories are searched recursively for `*.fastq`, `*.fastq.gz`, `*.fq`, and `*.fq.gz`. All chunks belonging to one amplicon are streamed into a single staged `.fastq.gz`; the originals are never modified.

## Current ABO configuration

Primer-derived target coordinates are stored in `config/samples.tsv`:

```text
sample  gene  amplicon   target_region                  phasing_min_length  core_depth_threshold
P001    ABO   fragment1  NG_006669.2:4228-18285        8000                50
P001    ABO   fragment2  NG_006669.2:11479-24671       8000                50
P001    ABO   fragment3  NG_006669.2:19001-32388       8000                20
```

The three overlapping products merge to the BED interval `NG_006669.2 4227 32388`.

## Important v0.6 outputs

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
```

The consensus FASTAs are reference-guided haplotype consensuses, not de novo assemblies. Sequence outside the primer-defined target, positions below the haplotype-specific callable-depth threshold, and unresolved variant spans are masked with `N`. An accepted deletion can exempt its deleted positions from the low-base-depth mask only on its ALT haplotype, with at least `callable_min_depth` exact ALT reads, a callable anchor, and no outside-target or uncertainty mask anywhere in its REF span. This permits supported deletions with zero aligned bases at deleted positions to be applied. For an unresolved insertion, the reference anchor is masked and the complete ambiguity remains documented in `uncertain_variants.tsv`.

Run focused regressions in an environment with Python, pysam, samtools, bcftools and tabix:

```bash
python -m unittest discover -s test -p 'test_*.py' -v
python test/integration_haplotype_consensus.py
```

## Running the workflow

From the repository root:

```bash
conda activate bloodgroup
snakemake -n -p
snakemake --cores 8 --software-deployment-method conda -p
```

Useful QC views:

```bash
column -t -s $'\t' results/qc/phasing/P001__ABO/P001__ABO.phasing_qc.tsv
column -t -s $'\t' results/qc/haplotypes/P001__ABO/P001__ABO.haplotag_qc.tsv
column -t -s $'\t' results/qc/haplotypes/P001__ABO/P001__ABO.coverage_qc.tsv
column -t -s $'\t' results/qc/haplotypes/P001__ABO/P001__ABO.variant_support.tsv | less -S
```

## Reference

The current ABO validation uses RefSeqGene `NG_006669.2` stored as `resources/references/ABO_reference.fasta`. Coordinates reported by this validation are RefSeqGene coordinates, not GRCh38 coordinates. Other genes and amplicons can be added as additional rows in `config/samples.tsv`.
