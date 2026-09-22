# KEL validation record — 2026-09-21

This record distinguishes executable checks from biological validation. The
repository clone is v0.6 at `ec722ae`. Local changes/configuration are uncommitted.

| Stage | Evidence/status |
| --- | --- |
| Repository inspection | README, Snakefile, all rule/env definitions, scripts, historical validation and smoke records inspected. |
| Native dependency resolution | PASS: core 202 packages and separate Clair3 134 packages, all `osx-arm64` or `noarch`; exact explicit package locks retained. |
| Core runtime | PASS: tool versions/imports in `tool_versions.json`. |
| Clair3 runtime/model format | PASS: native binaries; PyTorch HAC v5.2 checkpoints load and yield finite CPU outputs, `clair3_preflight.json`. |
| Original input integrity | PASS: 36 files, complete gzip CRC and FASTQ structure scan, SHA-256 inventory in `input_audit/`. |
| Reference and primers | PASS: NCBI `NG_007492.3`, 28,313 bp; all four primer pairs match exactly, `reference_provenance.json`. |
| Synthetic core workflow | PASS: 46 jobs completed; 48/48 expected primary reads mapped and passed MAPQ/length filters; expected BED matched (`smoke/verification.json`). |
| Synthetic Clair3 execution | PASS with limited scope: valid indexed empty VCF, 1,001 bp at 48x. No candidates; full-alignment calling not reached (`clair3_smoke/summary.json`). |
| Synthetic Clair3 candidates | PASS: 3/3 injected SNVs and expected genotypes recovered, including pileup and full-alignment inference (`clair3_candidate_smoke/summary.json`). |
| KEL core dry run | PASS: 46 planned jobs before analysis; `logs/kel-core-dry-run.log`. |
| KEL mapping/QC execution | PASS: all 46 jobs completed in 23:06; `logs/kel-core-run.log`. |
| KEL core artifact verification | PASS: BAM/index/filter/count/target/manifest checks and all 36 original SHA-256 hashes unchanged (`kel_core_verification.json`). |
| KEL calling/phasing dry run | PASS: 7 remaining jobs, `logs/kel-calling-dry-run.log`. |
| KEL calling/phasing execution | PASS: all 7 jobs completed in 3:29; `logs/kel-calling-run.log`. |
| KEL VCF verification | PASS: indexes, reference alleles, sample label, target bounds and record preservation (`kel_calling_verification.json`). |
| Calling/phasing final dry run | PASS: `logs/kel-final-dry-run.log` (prior stage). |
| Indel/phase/mask regression tests | PASS: 38 tests, `indel_fix/unit-tests.log`. |
| Synthetic support-to-consensus integration | PASS: exact known HP1/HP2 sequences and blocked deletion-mask case, `indel_fix/portable_synthetic/verification.json`. |
| KEL reconstruction dry run and execution | PASS: 11 jobs; final mask fix reran five affected jobs after another dry run, `indel_fix/kel-reconstruction-run.log` and `kel-mask-final-run.log`. |
| KEL reconstruction verification | PASS: read partitions, support accounting, masks, complete consensus sequence and 36 original hashes, `indel_fix/kel_reconstruction_verification.json`. |
| Leave-one-out tagging diagnostics | Both 560 CCT>C and 714 G>C remain ACCEPT when withheld from tagging; `indel_fix/leave_one_out/verification.json`. |
| Full workflow final dry run | PASS: nothing to do with all three KEL configs, `indel_fix/kel-final-dry-run.log`. |

The user confirmed one biological sample and requested no sample identifier.
The schema requires a grouping label, so both amplicons use `KEL`, yielding
analysis key `KEL__KEL`. This is a technical label only.

| Input | Barcode | FASTQ files | Reads | Bases | Median length | N50 |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Fragment 1 | 84 | 19 | 81,099 | 466,431,769 | 4,450 | 9,031 |
| Fragment 2 | 72 | 17 | 318,004 | 1,637,931,900 | 3,525 | 9,238 |

All 399,103 headers specify `dna_r10.4.1_e8.2_400bps_hac@v5.2.0`; the user also
confirmed this model. The native Clair3 environment includes converted PyTorch
`r1041_e82_400bps_hac_v520` checkpoints; checksums are recorded in the preflight.
No FASTQ filtering/downsampling is enabled. MAPQ 30 is used for derived variant
and phasing BAMs; only phasing adds the upstream 8 kb read-length threshold.
Before mapping, 21,998 (27.12%) fragment1 reads and 70,882 (22.29%) fragment2 reads
are at least 8 kb. Threshold suitability must be judged from actual mapping and
phase connectivity, not from the historical ABO thresholds alone.

The supplied fragment2 filenames lack chunk indices 2 and 5 within the observed
0–18 range. This does not demonstrate corruption or missing reads; export
completeness cannot be established without the original export/run manifest.
An additional full header scan found 399,103 unique read identifiers and no
duplicates across the supplied chunks (`input_audit/read_identifier_check.json`).

## Observed KEL core results

| Metric | Fragment 1 | Fragment 2 | Combined |
| --- | ---: | ---: | ---: |
| Input reads | 81,099 | 318,004 | 399,103 |
| Primary mapped reads | 74,571 | 279,035 | 353,606 (88.60%) |
| MAPQ >=30 primary mapped reads | 54,548 | 249,972 | 304,520 |
| Phasing reads (also length >=8 kb) | 15,484 | 70,406 | 85,890 |
| Minimum depth in designated target | 75 | 5,711 | 5,713 across union |
| Mean depth in designated target | 14,012.51 | 92,970.74 | 54,988.95 across union |
| Target bases below 50x | 0 | 0 | 0 |
| Fraction of depth outside designated amplicon | 0.628% | 0.913% | — |

Depth here comes from primary/MAPQ-filtered variant BAMs via the workflow's
`samtools depth -aa` outputs. Combined values sum both amplicons at each reference
position, including the overlap. These are read depths, not independent molecule
counts or haplotype-specific callability. Small fractions of depth extend outside
each designated product; no specific cause has been established. The supplied
primer targets, rather than observed depth intervals, define the calling BED.
The complete raw/filtered summaries are in `../workflow/results/summary/`.

## Observed KEL calling and phasing results

- Clair3 produced 46 raw records. Normalization/target filtering removed one
  record at `NG_007492.3:28038` outside the requested end coordinate 28023, leaving
  45 records. This is confirmed in `kel_target_filter_audit.json`.
- The normalized catalogue contains 44 biallelic records and one multiallelic
  `1/2` record. All 45 have the caller's PASS filter. PASS is a caller label, not
  independent confirmation of accuracy.
- The biallelic catalogue contains 41 heterozygous and 3 homozygous-alt records.
  All 41 eligible heterozygous records phased; none remained unphased.
- One observed phase set: `PS=516`, positions `516-25592`, inclusive span 25,077 bp.
  Thus this run's biallelic phase output connects the two amplicons. It does not
  establish switch-error rate without external truth or other validation.
- WhatsHap found 83,169 reads covering variants, kept 71,827 covering at least
  two variants, and selected 23 phase-informative reads for its phasing algorithm.
  The input phasing BAM contains all 85,890 reads that passed the workflow filters.
- The multiallelic site is retained in the normalized catalogue and excluded from
  the biallelic phasing input by design; it is not a phased call in this run.

Actual Clair3 inference processed 134 pileup candidates in 136.79 seconds and
132 full-alignment candidates in 27.26 seconds. The installed v2.0.3 C pileup
limit is 1,048,576 reads; the internal depth-144 parameter scales pileup tensor
channels when depth exceeds 216, rather than capping input to 144 reads. Observed
pileup candidate depths reached 148,118 (median 62,667). Full-alignment tensors
have their own model input dimensions. Thus the raw/filtered BAM depth should not
be confused with the model tensor's dimensions or haplotype-specific support.
See the [versioned Clair3 source](https://raw.githubusercontent.com/HKU-BAL/Clair3/v2.0.3/src/clair3_pileup.c)
and installed `envs/clair3-arm64/bin/clair3/CallVariantsFromCffi.py`.

Outputs:

```text
../workflow/results/variants/KEL__KEL/KEL__KEL.norm.vcf.gz
../workflow/results/variants/KEL__KEL/KEL__KEL.phasing_ready.vcf.gz
../workflow/results/phasing/KEL__KEL/KEL__KEL.phased.vcf.gz
../workflow/results/qc/phasing/KEL__KEL/KEL__KEL.phasing_qc.tsv
```

Calling/phasing execution is complete. After the requested indel repair,
haplotagged BAMs, support tables, masks and both consensus FASTAs were generated
and checked as described below. These are execution and internal consistency
checks, not a truth-set benchmark of the assay.

Targets use the user-requested conservative outer boundaries:

```text
Fragment 1: NG_007492.3:411-15302
Fragment 2: NG_007492.3:14079-28023
Union BED:  NG_007492.3  410  28023
Union: 27,613 bp; overlap: 1,224 bp
```

These are RefSeqGene coordinates, not GRCh38 coordinates. Reference SHA-256:
`488d6efe397c2fc5da5ecc65f43453cf6bc2e465d1c2e9beb8a55eaa5fc0074c`.
Alignment to this locus reference assesses consistency with KEL; it does not
independently measure specificity against the whole human genome.

## Historical evidence and limits

The repository documents ABO validation at v0.5/v0.6: 67 normalized records,
52/52 phased heterozygous sites in one block, 37,126 haplotagged reads, and four
unresolved candidates. These are historical narrative records in
`../workflow/docs/validation/`; their underlying BAM/VCF/consensus files and
execution logs are absent from the clone. They were not independently reproduced
here. The checked-in v0.3 smoke record explicitly says the DAG and tools were not
run. No previous KEL validation artifacts were present.

The new core smoke test verifies execution, mapping counts, filtering and target
generation against a deterministic artificial fixture. It is not a truth-set
assessment of variant sensitivity, genotype accuracy, phasing, or haplotypes.
An additional synthetic BAM copy with one homozygous-alt and two heterozygous
SNVs exercised Clair3 pileup, phasing, full-alignment and merging. All three exact
alleles/genotypes were recovered; this small engineered case is execution evidence,
not a representative variant benchmark. The source synthetic BAM stayed unchanged.
Clair3 2.0.3's PyTorch implementation differs from older TensorFlow versions;
model-name agreement alone does not establish equivalence to historical ABO calls.
The installed runner accepts the original amplicon flags but reports capping
`--ref_pct_full=1` to 0.3 internally; this behavior is retained in execution logs.

## Reconstruction defects repaired and revalidated

1. The original `workflow/scripts/haplotype_variant_support.py` counted deletion placeholders,
   reference skips, ambiguous bases and mismatching anchors as exact REF for
   indels. It also accepts matching indel events at an incorrect anchor as ALT.
   `indel_support_audit.json` reproduces four mismatching cases and two controls
   by extracting the original functions without running workflow side effects.
   The original script is archived in `indel_fix/before/`. The corrected counter
   requires an exact aligned allele window spanning the anchor, repeat context
   and right flank. Partial/conflicting observations are OTHER. SNV pileup
   handling also excludes placeholders and competing indel events from exact
   REF/ALT counts. All support thresholds are unchanged.
2. HP1/HP2 splitting previously used HP tags without requiring a single connected phase set.
   Whole-locus consensus could therefore join independently oriented phase blocks.
   Physical amplicon overlap alone does not establish haplotype connectivity.
   KEL did produce one phase set in this run, so disconnected blocks were not
   observed here. A new fail-fast guard now requires one informative phase set
   per target contig, with a recorded validation TSV.
3. The synthetic integration fixture exposed a second indel defect: accepted
   deletions produce no aligned bases at deleted positions, and the low-depth
   mask prevented bcftools from applying them. Deleted positions may now be
   exempted only on the accepted ALT haplotype, with at least 50 exact ALT reads,
   a callable anchor, and no uncertainty/outside-target mask anywhere across the
   REF allele. If any part blocks the deletion, its low-depth positions stay
   masked. Both successful deletion application and blocked-anchor behavior have
   exact-sequence integration tests.

`config/kel.reconstruction.yaml`, applied after the core and calling configs,
enables the repaired stages. All 38 unit tests and the portable synthetic
integration test passed. The reconstruction dry run planned 11 jobs; all ran
successfully in 53.8 seconds. After the final mask precedence fix, the five
affected final jobs were dry-run and completed again.

| Result | HP1 | HP2 |
| --- | ---: | ---: |
| Assigned reads | 32,100 | 21,523 |
| Minimum target depth | 630 | 13 |
| Mean target depth | 13,371.93 | 8,762.21 |
| Target positions below 50x | 0 | 87 (411–497) |
| Final FASTA length | 28,313 | 28,311 |
| N bases | 724 | 811 |
| Masked target reference positions | 24 | 111 |

There are 32,267 unassigned reads; assigned plus unassigned exactly account for
all 85,890 phasing reads. All assigned reads have PS516. The support table has
33 accepted variants (32 SNVs and `560 CCT>C`) and 12 unresolved variants at
624, 639, 660, 692, 696, 717, 727, 758, 759, 1688, 25592 and 27258. All unresolved
REF spans are masked in both FASTAs (24 reference bases). The 700 bases outside
the target are masked in both; HP2 additionally has 87 low-depth target bases.
The accepted two-base deletion is applied in HP2. No KEL deletion required the
zero-depth mask exemption; that behavior is tested by the synthetic fixture.

The archived defective counter was rerun in a separate scratch directory on
identical KEL HP BAMs. One decision changes: `759 CA>C`, previously ACCEPT,
becomes UNRESOLVED. Its HP2 REF/ALT/OTHER counts change from `99/860/1` to
`4/4/952` when the complete allele window is required. This is insufficient exact
support, not proof the biological deletion is absent. The accepted `560 CCT>C`
retains 757 exact HP2 ALT reads, ALT fraction 0.8885 and OTHER fraction 0.1125.
All 44 biallelic records have count changes, including the stricter SNV handling;
the other status decisions agree. See `indel_fix/kel_legacy_comparison.json`.

An independent verifier checked the full read partition, every catalogue/support
record, accepted VCF alleles/genotypes, indel denominators, 200 aligned-pair read
windows, masks and complete consensus sequences. Both sequences exactly match
independent reference-coordinate reconstruction. All 36 raw SHA-256 hashes were
rechecked and remain unchanged. See `indel_fix/kel_reconstruction_verification.json`
and the versioned report `../workflow/docs/validation/KEL_indel_revalidation.md`.
This completes computational repair and revalidation of the available KEL run;
external variant/phase truth and representative assay validation remain absent.

Leave-one-out checks separately withheld accepted deletion `560 CCT>C` and the
weakest accepted SNV by expected ALT fraction, `714 G>C`, from haplotagging while
retaining original phased genotypes for support evaluation. Both remain ACCEPT:
HP2 ALT fractions 0.7218 and 0.3841; between-HP differences 0.7107 and 0.3794.
All alignment identities and PS516 were preserved. Withholding 560 changed
scratch `727 G>A` from UNRESOLVED to ACCEPT, demonstrating assignment sensitivity;
its production mask remains unchanged. Withholding 714 changed no status.
Production VCF/support hashes are unchanged. These checks retain the original
phase solution, so they are not independent phase truth. See
`indel_fix/leave_one_out/README.md` and `verification.json`.

## Local workflow changes

- Mark `all` as the default target and add explicit `qc_only` for new-locus QC.
- Allow `qc.nanoplot_no_static` to produce statistics/HTML without browser exports.
- Honor an explicit Clair3 model path and quote its executable, enabling a
  separate native caller environment.
- Add local KEL sample/config files; retain the default ABO configuration.
- Repair exact indel/SNV support and accepted-deletion masking; add a phase-set
  guard, `pysam` haplotype dependency, regression tests and KEL validation report.

See `README.md` for runnable commands and environment details. Generated outputs
remain under `../workflow/results/`; original sequencing files remain outside
the clone and are accessed read-only.
