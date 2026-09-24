# P001 ABO validation — v0.6 haplotagging and haplotype reconstruction

This document records the real-data validation used to design the v0.6 haplotype reconstruction stage. It extends the v0.5 phasing validation rather than replacing it.

## Starting point

- Analysis: `P001__ABO`
- Reference: RefSeqGene `NG_006669.2`
- Target union: `NG_006669.2:4228-32388` (28,161 bp)
- Three independently sequenced overlapping long-range PCR amplicons
- Complete normalized catalogue: 67 variant records
- Phasing-ready biallelic catalogue: 66 records
- WhatsHap-phased heterozygous variants: 52/52
- Phase sets: one (`PS=9033`)
- Phase block: positions 9033-29201 (20,169 bp inclusive)

## Haplotagging validation

The validated phased VCF and long-read-enriched phasing BAM were passed to `whatshap haplotag`.

Read assignment:

```text
Input phasing BAM reads    42,373
HP1 reads                  18,412
HP2 reads                  18,714
Assigned reads             37,126
Unassigned reads            5,247
Assigned fraction           87.6%
HP1 among assigned          49.6%
HP2 among assigned          50.4%
```

All 37,126 haplotagged reads carried `PS=9033`, matching the single validated phase block.

## Haplotype-specific target coverage

Coverage was measured separately after splitting the haplotagged BAM by the `HP` tag.

```text
                         HP1          HP2
Target positions       28,161       28,161
Mean depth            7,715.68x    7,727.58x
Minimum depth            71x          527x
Zero-depth positions       0            0
Positions <10x             0            0
Positions <20x             0            0
Positions <50x             0            0
Positions <100x            4            0
```

At the v0.6 validation threshold of 50x, both haplotypes were callable across 100% of the primer-defined target.

## Multiallelic repeat at position 9109

The complete Clair3 catalogue contains:

```text
POS   REF               ALT    QUAL  GT
9109  CATATATATATATAT   C,CAT  7.26  1/2
```

WhatsHap does not routinely phase this multiallelic `1/2` record. CIGAR-level deletion counting in independently haplotagged reads nevertheless showed haplotype-associated repeat-length distributions.

Most common deletion events:

```text
HP1
12 bp deletion   4,924
14 bp deletion   3,273
10 bp deletion   1,483
16 bp deletion     858

HP2
14 bp deletion   5,818
16 bp deletion   2,097
12 bp deletion   1,948
10 bp deletion     380
```

The dominant evidence is consistent with HP1 carrying the 12-bp deletion allele and HP2 the 14-bp deletion allele, i.e. a probable `2|1` orientation relative to the Clair3 ALT order. However, the broad repeat-length distributions make the exact allele insufficiently clean for automatic consensus insertion. The site is therefore retained as `UNRESOLVED` and its REF span is masked.

## Haplotype-specific support audit

All 66 biallelic records were evaluated using exact REF/ALT observations in HP1 and HP2 reads. Most SNVs and indels showed strong reciprocal support consistent with the WhatsHap phase.

Examples of strongly supported homozygous-alt SNVs included:

```text
5147   T>C   HP1 ALT 0.9918   HP2 ALT 0.9929
5165   G>A   HP1 ALT 0.9740   HP2 ALT 0.9729
5801   T>C   HP1 ALT 0.9881   HP2 ALT 0.9869
27120  G>A   HP1 ALT 0.9895   HP2 ALT 0.9827
28096  A>G   HP1 ALT 0.9871   HP2 ALT 0.9874
29745  T>C   HP1 ALT 0.9409   HP2 ALT 0.9433
```

Several low-QUAL heterozygous variants remained clearly haplotype-associated. For example, `28256 T>C` had Clair3 QUAL 6.43 but HP1 ALT fraction 0.1621 and HP2 ALT fraction 0.8175. This demonstrated that a hard global QUAL threshold would remove genuine haplotype-associated signal.

## Weak heterozygous sites and leave-one-out validation

Three phased heterozygous SNVs had an HP1/HP2 ALT-fraction difference below 0.40:

```text
POS    PHASE   HP1 ALT   HP2 ALT   delta
10219  1|0     0.3814    0.0177    0.3637
19755  1|0     0.6476    0.3478    0.2998
20067  0|1     0.3356    0.6630    0.3274
```

Because these variants themselves contributed to the original haplotagging, each site was removed from the phased VCF and haplotagging was repeated. The tested position was then measured in reads assigned using the remaining phased variants.

Leave-one-out results:

```text
POS    PHASE   HP1 ALT   HP2 ALT   delta
10219  1|0     0.3721    0.0301    0.3420
19755  1|0     0.6343    0.3611    0.2732
20067  0|1     0.3476    0.6504    0.3028
```

All three retained the direction predicted by their original phase after their own contribution to haplotagging was removed. They were therefore accepted for reconstruction. This empirical result motivated a configurable minimum haplotype ALT-fraction difference of 0.25 rather than an arbitrary high QUAL cutoff.

## Unresolved biallelic indels

Three nominal homozygous-alt indels conflicted strongly with haplotype-specific exact allele support:

```text
POS    REF>ALT       GT   QUAL   HP1 ALT   HP2 ALT   interpretation
13629  CA>C          1/1   7.13   0.5032    0.4958    competing repeat-associated indels
13963  CT>C          1/1   0      0.4561    0.4604    LowQual + competing indels
25292  G>GACATACAC   1/1  11.58   0.0227    0.1579    exact ALT poorly supported
```

At 13629 and 13963, the number of `OTHER` indel observations exceeded the exact REF or ALT counts on both haplotypes. At 25292, most observations were competing indels and the claimed exact insertion had little support. These sites are classified as `UNRESOLVED`, not as confidently reference or confidently alternate.

## v0.6 classification strategy

The validation supports the following conservative reconstruction rules:

- preserve the complete multiallelic normalized VCF as the variant catalogue;
- use the biallelic phased VCF for HP1/HP2 read assignment;
- require adequate haplotype-specific read depth;
- accept phased heterozygous variants when the expected ALT haplotype has sufficient exact ALT support and exceeds the opposite haplotype by a configurable margin;
- accept homozygous-alt variants only when both haplotypes independently support the exact ALT allele;
- do not use a hard global QUAL threshold;
- keep non-PASS, multiallelic, unsupported complex, unphased, or support-conflicting candidates as `UNRESOLVED`;
- mask unresolved spans, low-depth positions and sequence outside the primer-defined target with `N` in the reference-guided haplotype FASTAs.

The v0.6 defaults derived from this validation are `callable_min_depth=50`, `min_support_depth=20`, `min_het_alt_fraction=0.30`, `min_het_delta=0.25`, `min_hom_alt_fraction=0.80`, and `max_other_fraction=0.25`. These values remain configuration parameters and should be revalidated when the assay, chemistry, basecaller, amplicon design, or application changes materially.


## Consensus round-trip validation

The final reference-guided HP1 and HP2 FASTAs were independently remapped to `NG_006669.2` with minimap2 using an assembly-to-reference alignment preset. Each haplotype produced exactly one primary alignment, with no secondary or supplementary alignments.

Expected variants were derived directly from the consensus-ready phased VCF and compared with variants recovered from the remapped haplotype FASTAs.

```text
HP1 expected variants      48
HP1 SNVs recovered         45/45
HP1 indels confirmed        3/3
HP1 unexpected variants     0

HP2 expected variants      26
HP2 SNVs recovered         24/24
HP2 indels confirmed        2/2
HP2 unexpected variants     0
```

The five indels not re-called by haploid `bcftools mpileup/call` from the single synthetic consensus alignment were confirmed directly in the alignment/pileup at the expected positions and with the expected sequence changes:

```text
HP1
24742  CACAG>C
25504  CACAG>C
28334  T>TGAGGC

HP2
15704  A>ACAGTTTGG
24543  AC>A
```

Thus, all expected accepted alleles were present on the intended reconstructed haplotype and no unexpected sequence variants were introduced by the consensus-generation step.

## Biological interpretation and sample origin

The phased ABO coding pattern contains a coherent A-like haplotype and a coherent B-like haplotype. The B-associated coding markers occur with the same phase orientation, supporting a conventional AB-type ABO genotype for P001.

P001 originates from the Promega Human Genomic DNA control used for this validation dataset.

This biological interpretation is recorded separately from the computational validation: the workflow reconstruction remains based on the sequence evidence, while the Promega control provides the sample provenance for the P001 validation run.
