# Long-range ONT Snakemake pipeline — v0.2

This version matches the current ABO dataset structure: **one patient + one gene + multiple long-range PCR amplicons**, where each ONT run folder contains many FASTQ chunks.

## What changed from v0.1

- Each row now represents one **amplicon** (e.g. ABO fragment 1/2/3).
- `fastq_input` may be a **directory or a single FASTQ file**.
- Directories are searched recursively for `.fastq`, `.fastq.gz`, `.fq`, and `.fq.gz`.
- Compressed and uncompressed files can be mixed.
- All chunks for an amplicon are streamed into one gzipped staging FASTQ; the original files are never modified.
- QC and mapping are retained per amplicon.
- Amplicon BAMs are merged at the **sample + gene** level before variant calling/phasing.
- An `input_manifest.tsv` records every original FASTQ chunk used in the analysis.

## Workflow

```text
FASTQ chunks (.fastq and/or .fastq.gz)
        │
        ├─ ABO fragment 1 ─ combine ─ NanoPlot ─ minimap2 ─ BAM ─┐
        ├─ ABO fragment 2 ─ combine ─ NanoPlot ─ minimap2 ─ BAM ─┼─ samtools merge
        └─ ABO fragment 3 ─ combine ─ NanoPlot ─ minimap2 ─ BAM ─┘
                                                                  │
                                                          merged ABO BAM
                                                                  │
                                                   mapping + coverage QC
                                                                  │
                                            [optional] Clair3 → WhatsHap
                                                                  │
                                                 haplotype consensus
```

## Configuration for your present ABO data

Edit `config/samples.tsv`:

```tsv
sample  gene  amplicon   fastq_input                                      reference
P001    ABO   fragment1  data/ABO_Fragment_1-20260911T075013Z-1-001      resources/references/ABO_reference.fasta
P001    ABO   fragment2  data/ABO_Fragment_2-20260911T075014Z-1-001      resources/references/ABO_reference.fasta
P001    ABO   fragment3  data/ABO_Fragment_3-20260911T075016Z-1-001      resources/references/ABO_reference.fasta
```

Tabs are required in this file.

You may also use absolute paths. On WSL, a Windows path such as `C:\Users\...\ABO_Fragment_1...` is normally addressed as `/mnt/c/Users/.../ABO_Fragment_1...`.

## Reference

For the first real test, place one reference FASTA at:

```text
resources/references/ABO_reference.fasta
```

All ABO amplicons in the same analysis must use the same reference. Prefer a genomic ABO target/locus reference containing the regions covered by all three PCR fragments.

## Recommended Windows setup

The workflow tools (Snakemake, minimap2, samtools, NanoPlot, Clair3, WhatsHap) are easiest to run under **WSL2/Linux** with Miniforge/Mambaforge rather than directly in Windows PowerShell.

Example:

```bash
mamba create -n bloodgroup-smk -c conda-forge -c bioconda snakemake
conda activate bloodgroup-smk
```

From the pipeline directory:

```bash
snakemake -n -p
snakemake --use-conda --cores 8 -p
```

## First real test: QC + mapping only

The provided `config/config.yaml` has variant calling/phasing/consensus **disabled by default**. This is intentional. First verify the three amplicons map as expected.

Expected key outputs:

```text
results/summary/input_manifest.tsv
results/summary/amplicon_summary.tsv
results/summary/gene_summary.tsv
results/qc/raw/P001__ABO__fragment1/NanoStats.txt
results/mapping/amplicons/P001__ABO__fragment1/...
results/mapping/genes/P001__ABO/P001__ABO.merged.bam
results/mapping/genes/P001__ABO/P001__ABO.coverage.txt
```

Once mapping/QC looks correct, configure the correct Clair3 ONT model and enable:

```yaml
workflow:
  run_variant_calling: true
  run_phasing: true
  run_consensus: true
```

## Synthetic smoke test included

The `test/` directory contains three synthetic amplicon folders using both uncompressed and gzipped FASTQs. Its purpose is to test the new universal input/combination layer. The full mapping workflow still requires Snakemake/minimap2/samtools/NanoPlot to be installed.
