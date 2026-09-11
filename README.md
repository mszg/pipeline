# Long-range ONT Snakemake pipeline — v0.3

This workflow is for targeted Oxford Nanopore long-range amplicon sequencing. It currently supports the real ABO test dataset used during development: one patient, one gene, three overlapping PCR amplicons, and many FASTQ chunks per amplicon.

## v0.3 patch

v0.3 keeps the raw alignment outputs from v0.2, but adds an explicit separation between **raw mapping**, **small-variant analysis**, and **phasing**.

The change was motivated by the first real ABO run. Fragment 1 contained many short/truncated reads and large numbers of secondary/supplementary alignments, while reads approaching the expected full amplicon length mapped much more confidently. The pipeline now measures that behavior instead of treating every alignment record as equivalent.

### New in v0.3

- reports **primary mapped reads** rather than relying only on raw `samtools flagstat` mapping percentages;
- reports primary mappings at **MAPQ >=20, >=30 and >=50**;
- reports **secondary and supplementary alignment counts**;
- produces a **read-length vs MAPQ** QC table for every amplicon;
- creates a dedicated **variant BAM** containing only primary, mapped, MAPQ-filtered alignments;
- creates a separate **phasing BAM** that additionally enriches for long reads;
- supports an optional per-amplicon `target_region` such as `NG_006669.2:4000-18500`;
- calculates contiguous **high-depth core intervals** from the cleaned variant BAM;
- keeps all raw BAMs intact for traceability and future structural-variant analyses;
- Clair3 now consumes the cleaned gene-level variant BAM;
- WhatsHap now consumes the long-read-enriched gene-level phasing BAM;
- the new BAM QC/filtering Python environment is pinned to Python 3.12 for reproducibility without changing the already validated minimap2/samtools environment.

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
     │
     ▼
  minimap2
     │
     ▼
RAW sorted BAM  ─────────────── retained unchanged
     │
     ├─ alignment QC
     │    ├─ primary mapped %
     │    ├─ MAPQ >=20/30/50
     │    ├─ secondary/supplementary counts
     │    └─ length-vs-MAPQ table
     │
     ├──────────────────────────────────┐
     ▼                                  ▼
VARIANT BAM                         PHASING BAM
primary only                        primary only
mapped only                         mapped only
MAPQ >=30                           MAPQ >=30
no length cutoff                    long-read enriched (default >=8 kb)
optional target region              optional target region
     │                                  │
     ├─ coverage/depth QC               │
     ├─ core interval detection         │
     │                                  │
     ▼                                  ▼
merge by patient+gene              merge by patient+gene
     │                                  │
     ▼                                  ▼
Clair3                              WhatsHap
     │                                  ▲
     └──────────── VCF ─────────────────┘
```

Clair3/WhatsHap remain disabled by default until the correct ONT Clair3 model is configured.

## Current ABO configuration

`config/samples.tsv` contains:

```tsv
sample  gene  amplicon   fastq_input                                      reference                               target_region  phasing_min_length  core_depth_threshold
P001    ABO   fragment1  data/ABO_Fragment_1-20260911T075013Z-1-001       resources/references/ABO_reference.fasta                 8000                50
P001    ABO   fragment2  data/ABO_Fragment_2-20260911T075014Z-1-001       resources/references/ABO_reference.fasta                 8000                50
P001    ABO   fragment3  data/ABO_Fragment_3-20260911T075016Z-1-001       resources/references/ABO_reference.fasta                 8000                20
```

The file is tab-separated. `target_region` is deliberately blank for the present dataset. The intervals observed from coverage should **not** automatically be treated as the expected PCR coordinates. Once the exact primer-derived amplicon coordinates are available, they can be entered as, for example:

```text
NG_006669.2:4225-18287
```

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

It deliberately does **not** apply a minimum read-length cutoff by default, because shorter reads can still contain valid local SNP/indel evidence.

### Phasing BAM

The phasing BAM uses the same alignment cleanup and additionally requires a configurable minimum read length. This is intended to enrich for reads carrying long-range linkage information. It does not replace the raw BAM.

## Important output files

```text
results/summary/input_manifest.tsv
results/summary/amplicon_summary.tsv
results/summary/gene_summary.tsv

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

Do not enable Clair3 until `clair3.model_path` points to the correct ONT model for the basecalling chemistry/model used for the sequencing run.

Once configured:

```yaml
workflow:
  run_variant_calling: true
  run_phasing: true
  run_consensus: true
```

Clair3 will use:

```text
results/mapping/genes/P001__ABO/P001__ABO.variant.bam
```

and WhatsHap will use:

```text
results/mapping/genes/P001__ABO/P001__ABO.phasing.bam
```

## Reference

The current ABO test uses the full RefSeqGene sequence `NG_006669.2` saved as:

```text
resources/references/ABO_reference.fasta
```

The workflow itself is not ABO-specific. Other genes/amplicons can be added as additional rows in `config/samples.tsv`.
