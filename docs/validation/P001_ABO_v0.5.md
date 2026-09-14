# P001 ABO validation — v0.5 phasing

This document records the real-data validation used to design the v0.5 phasing stage of the long-range ONT amplicon pipeline.

## Dataset

- Analysis: `P001__ABO`
- Reference: RefSeqGene `NG_006669.2`
- Three independently sequenced, overlapping long-range PCR amplicons:
  - fragment1: `NG_006669.2:4228-18285` (14,058 bp)
  - fragment2: `NG_006669.2:11479-24671` (13,193 bp)
  - fragment3: `NG_006669.2:19001-32388` (13,388 bp)
- Merged target BED: `NG_006669.2 4227 32388`
- ONT basecalling model: `dna_r10.4.1_e8.2_400bps_hac@v5.2.0`
- Clair3 model: `r1041_e82_400bps_hac_v520`

## Independent variant-signal check

Before accepting Clair3 output, the cleaned variant BAM was tested independently with `bcftools mpileup/call`. This produced 69 candidate records, including multiple high-depth, approximately 50:50 heterozygous SNPs and near-homozygous alternate SNPs. This demonstrated that clear non-reference signal was present in the BAM and helped isolate the initial zero-call problem to Clair3 configuration rather than mapping or target definition.

Removing `--chunk_num=-1` from the Clair3 invocation resolved the zero-call behavior.

## Clair3 variant calling

The target-aware Clair3 run produced 68 records before target-restricted unsplit normalization. The previous v0.4 normalization used `bcftools norm -m -any`, which split a multiallelic repeat-associated site at position 9109 into duplicate genomic positions.

Original Clair3 representation at 9109:

```text
POS   REF               ALT    QUAL  GT   DP     AD
9109  CATATATATATATAT   C,CAT  7.26  1/2  25598  1132,9782,7383
```

After split normalization, this became two biallelic records at the same coordinate. WhatsHap warned about the duplicate position and skipped one representation.

## Multiallelic-preserving validation

The VCF was therefore normalized without `-m -any`, retaining the multiallelic genotype structure. After target restriction, the complete normalized catalogue contained:

- 67 variant records/sites
- 66 biallelic records
- 1 multiallelic record (`NG_006669.2:9109`, genotype `1/2`, QUAL 7.26)

WhatsHap 2.8 did not phase the unsplit `1/2` multiallelic site. It remained unphased (`PS=.`). This site is in a repetitive AT tract and is preserved in the complete variant catalogue rather than being silently discarded or forced into an incorrect biallelic phase representation.

## Biallelic phasing-ready validation

A separate phasing-ready VCF was generated with:

```bash
bcftools view -m2 -M2
```

This produced 66 biallelic records with no duplicate genomic positions.

WhatsHap 2.8 results:

- usable heterozygous variants: 52
- phased heterozygous variants: 52/52 (100%)
- homozygous-alt biallelic variants: 14
- phase sets: 1
- phase-set ID: `9033`
- largest phase block: 52 variants
- phase-block interval: `NG_006669.2:9033-29201`
- reads covering the 52 variants: 37,820
- reads covering at least two variants: 37,199
- most phase-informative reads selected by WhatsHap: 38
- elapsed time: 11.9 s
- maximum memory: 0.124 GB

Genotype summary of the phased biallelic VCF:

```text
15  0|1  PS=9033
37  1|0  PS=9033
14  1/1  PS=.
```

The `0|1` versus `1|0` orientation is arbitrary; the important result is that all 52 accessible heterozygous variants share one consistent phase set.

## v0.5 design decision

The validated architecture is:

```text
Clair3 raw VCF
    |
    v
multiallelic-preserving, target-restricted normalized VCF
    |\
    | \__ complete variant catalogue for reporting
    |
    v
biallelic phasing-ready VCF
    |
    v
WhatsHap
    |
    v
phased VCF + automated phasing QC
```

Multiallelic variants are retained in the complete catalogue but are not silently split into duplicate positions for routine WhatsHap phasing. Haplotype consensus/reconstruction remains disabled in v0.5 pending explicit handling and validation of unphased/multiallelic and low-confidence variants.
