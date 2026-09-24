# KEL Apple Silicon setup, indel repair and validation

The original indel counter could count ambiguous/deleted/skipped anchors as REF and accept an indel with the wrong anchor as ALT. The repair checks complete aligned allele windows, prevents deletion masks from blocking supported alleles, and rejects disconnected phase sets before reconstruction. KEL now completes through both masked haplotype FASTAs on native Apple Silicon, with executable validation evidence and reproducible setup assets.

For example, the old counter accepted `NG_007492.3:759 CA>C` with HP2 REF/ALT/OTHER counts `99/860/1`. Requiring a complete allele window changes those counts to `4/4/952`; the site is now unresolved and masked. The supported `560 CCT>C` deletion remains accepted and is applied to HP2. Acceptance thresholds were not relaxed.

## Workflow and reconstruction changes

### Exact allele support

- Add an importable, testable `allele_support.py` helper and an explicit `pysam` dependency in the haplotype environment.
- Count normalized pure insertions/deletions from each read's CIGAR and aligned query sequence across the VCF anchor, complete allele context and an aligned right flank.
- Extend the window through equivalent right-shifted placements in repeats so an alternate event represented elsewhere in the repeat cannot be mistaken for REF.
- Require exact window equality for REF or ALT. Wrong anchors, ambiguous bases, partial/clipped windows, missing flanks, reference skips and conflicting sequence count as OTHER.
- Exclude unmapped, secondary, supplementary, QC-failed and duplicate alignments from indel support. Count all remaining overlapping reads without a depth cap; the existing MAPQ/read-length filters still define the HP BAMs.
- Validate the VCF REF against the reference. Unsupported complex alleles and indels without usable flanking sequence remain unresolved. Nearby variation in the window can conservatively reduce exact support.
- Tighten the SNV pileup parser: handle start/end markers and multi-digit indel lengths; reject malformed input; classify deletion placeholders, skips, ambiguous bases and competing indel events as OTHER. Ignore zero-depth placeholder lines and avoid running mpileup on an empty SNV site list.
- Append `SUPPORT_METHOD` to support and uncertainty TSVs (`aligned_allele_window`, `pileup_base`, or `unsupported`). Track the helper script as a Snakemake input so its changes invalidate dependent support outputs.
- Preserve full normalized-catalogue accounting, accepted-variant VCF generation, uncertainty reasons and BED masks. Multiallelic candidates remain in the catalogue and are masked when unsupported by the reconstruction stage.

### Phase connectivity guard

- Add `validate_haplotype_phase_sets` as a prerequisite of haplotagging and HP splitting.
- Require a nonempty target BED, exactly one VCF sample matching the analysis key, and exactly one informative phase set per target contig.
- Fail on missing PS values, multiple disconnected blocks, a target contig without informative phased heterozygotes, unexpected contigs or malformed single-sample VCF input.
- Record phase set, phased/unphased heterozygote counts, start/end coordinates and PASS in `phase_set_validation.tsv`.
- Allow unphased heterozygotes to proceed to uncertainty masking. Separate contigs have independent HP orientations; this check does not establish chromosome-wide phase across contigs.

### Deletion-aware consensus masking

- Pass the support table, haplotype number and callable-depth threshold into each consensus-mask rule.
- Exempt deleted reference positions from the low-base-depth mask only when that deletion is accepted on the given ALT haplotype and exact ALT read support reaches the callable threshold.
- Require a callable anchor and no outside-target or uncertainty mask anywhere across the deletion's REF span. If an immutable mask or low-depth anchor blocks the deletion, retain its low-depth positions so skipped deletions cannot expose unsupported reference bases.
- Keep reference-haplotype, low-support and unresolved deletion positions masked. Outside-target and uncertainty masks always retain precedence.
- Continue using bcftools to produce reference-guided HP1/HP2 FASTAs. The synthetic integration test checks that accepted insertions/deletions actually change the sequence, and that blocked deletions retain the expected Ns.

### Workflow entry points and caller/QC integration

- Mark `all` explicitly as the default target.
- Add a reference-independent `qc_only` entry point for input staging, the input manifest and raw NanoPlot QC before target/reference setup.
- Add optional `qc.nanoplot_no_static`; enabled for KEL to produce statistics/HTML without browser-driven static plot export. Existing default behavior is retained when the option is absent.
- Add explicit `clair3.model_path` support, validate that the model directory exists, quote the configured executable, and retain executable-relative model discovery as a fallback.

## KEL configuration and native Apple Silicon setup

Three ordered configuration overlays allow each stage to be dry-run and validated separately:

| Configuration | Enabled behavior |
| --- | --- |
| `config/kel.yaml` | Input staging, raw QC, mapping, filtering, coverage, targets and summaries |
| `config/kel.calling.yaml` | Adds target-aware Clair3 calling and WhatsHap phasing through a separate native caller wrapper |
| `config/kel.reconstruction.yaml` | Adds guarded haplotagging, HP support/coverage, uncertainty masks and both consensuses |

- `kel.samples.tsv` groups two amplicons from one biological sample under the technical label `KEL`, producing `KEL__KEL`; no patient identifier was assigned.
- Fragment 1 uses barcode84 and `NG_007492.3:411-15302`; fragment 2 uses barcode72 and `NG_007492.3:14079-28023`. Union BED is `NG_007492.3 410 28023`: 27,613 bp with 1,224 bp overlap.
- Reference provenance records the NCBI RefSeqGene accession/version, 28,313 bp length, source URL, SHA-256 and exact checks for all four supplied primer pairs. Coordinates are RefSeqGene, not GRCh38.
- Match the observed/user-confirmed basecaller `dna_r10.4.1_e8.2_400bps_hac@v5.2.0` to the bundled `r1041_e82_400bps_hac_v520` model.
- Retain `map-ont`, MAPQ 30 for variant/phasing BAMs, an 8 kb minimum read length only for phasing, and two requested cores. KEL FASTQ filtering/downsampling stays disabled; original reads and raw mapped BAMs are preserved.
- Keep the original support thresholds: callable depth 50; exact support depth 20 per HP; expected heterozygous ALT fraction 0.30; between-HP ALT-fraction difference 0.25; homozygous ALT fraction 0.80; OTHER fraction at most 0.25; caller PASS required. The 100,000 mpileup limit applies to SNVs only.
- Preserve the existing default ABO configuration and historical validation records.

`deployment/kel-apple-silicon/` publishes the reusable setup assets previously stored beside the clone:

- Native compatibility assessment performed before installation, with source references and Python/architecture/model constraints.
- Separate core Python 3.12 and Clair3 Python 3.11 Conda YAML specifications, plus exact `osx-arm64`/`noarch` package locks. No existing Conda environment or Homebrew packages were modified during setup.
- Tested core versions: Snakemake 9.27.0, minimap2 2.31, samtools/bcftools/htslib 1.24, pysam 0.24.1, WhatsHap 2.8, NanoPlot 1.48.0 and filtlong 0.3.1. Native Clair3 2.0.3 uses bundled PyTorch models.
- A Snakemake wrapper that uses workspace caches, two requested cores, the greedy scheduler, incomplete-job reruns and one BLAS/OpenMP thread per process. A separate Clair3 wrapper selects its own PATH/CONDA_PREFIX and CPU runtime.
- A model preflight covering native binaries, checkpoint checksums, installed-loader compatibility, finite CPU inference for both models and parsing of the workflow's caller arguments.
- A read-only full gzip/FASTQ audit with original-byte SHA-256 hashes, structure/CRC validation, length/quality summaries, aggregate header metadata, chunk-index gap reporting and source size/mtime preservation.
- Real-output verifiers for core BAM/index/filter/target accounting, VCF/reference/sample/catalogue checks, HP read partitions, indel denominators, masks, complete FASTA reconstruction and unchanged originals.
- A deterministic 48-read core fixture; an empty-output Clair3 execution smoke; a three-injected-SNV Clair3 candidate/inference smoke; archived defective scripts for before/after comparison; and two reproducible leave-one-out tagging diagnostics.
- Selected original run logs and numerical JSON evidence. `MANIFEST.json` records source/published checksums and edits that replace historical machine paths or derive script roots from their deployed location. The documented deployment preserves the tested sibling `setup/` and `workflow/` layout.
- Original sequencing files, full input audits, installed environments, model binaries and generated BAM/VCF/FASTA datasets remain local under the existing data/output ignore policy.

## Pipeline features retained and exercised

The existing pipeline supports multi-file FASTQ discovery (including compressed chunks and input paths containing spaces), deterministic input staging and manifests, optional read filtering, NanoPlot QC, minimap2 long-read alignment, reference indexing, per-amplicon and merged gene BAMs, distinct variant/phasing filters, indexed BAM/VCF dependencies, primer-target BED generation and merging, alignment/length/MAPQ/depth QC, amplicon/gene summaries, target-aware Clair3 calling, normalization/catalogue retention, biallelic WhatsHap phasing, haplotagging, haplotype coverage/support reporting, uncertainty masks and reference-guided consensuses. This change wires KEL into that workflow and repairs the support/mask/connectivity behavior described above; it does not replace those established stages.

## Executed validation

| Check | Observed evidence |
| --- | --- |
| Original inputs | 36 valid gzip/FASTQ files, 399,103 reads; all source hashes unchanged; complete read-ID audit found no duplicates |
| Synthetic core | All 46 jobs completed; 48/48 expected reads mapped and passed configured filters; target BED matched |
| Synthetic Clair3 | Indexed empty VCF execution passed; separate candidate smoke recovered 3/3 injected SNVs with expected genotypes through pileup and full-alignment inference |
| KEL core | Dry run passed; all 46 jobs completed; 353,606 primary mapped reads and 304,520 variant-filtered reads |
| KEL calling/phasing | Dry run passed; all seven jobs completed; 46 raw records, one out-of-target record removed, 45 retained catalogue records |
| Phase connectivity | 41/41 eligible heterozygotes phased in one set, PS516, positions 516–25592; three homozygous ALT and one unphased multiallelic catalogue record retained |
| KEL reconstruction | Dry run planned 11 jobs; all completed in 53.8 seconds; final mask refinement dry-ran and reran five affected jobs |
| HP read accounting | 32,100 HP1 + 21,523 HP2 + 32,267 unassigned = all 85,890 input phasing reads; assigned reads retain PS516 |
| Support decisions | 33 ACCEPT (32 SNVs, one deletion); 12 UNRESOLVED; unchanged thresholds |
| HP1 FASTA | 28,313 bases, 724 Ns, including 24 masked target reference positions |
| HP2 FASTA | 28,311 bases, 811 Ns, including 111 masked target reference positions; accepted two-base deletion applied |
| Mask accounting | 700 outside-target positions in both; 24 uncertainty positions in both; HP2 additionally has 87 positions below 50x at 411–497 |
| Independent reconstruction | Both complete FASTA sequences exactly match independently applied alleles/masks; all indel denominators and 200 independently extracted read windows checked |
| Regression tests | All 38 unit tests pass: exact/shifted indels, partial windows, wrong anchors, skips/placeholders, long events, flags, malformed pileup, phase guards and mask precedence |
| Portable integration | Known SNV/insertion/deletion, wrong-anchor event, multiallelic/unphased uncertainty, low depth and outside-target masks; both complete synthetic sequences match, including a separate blocked-deletion-anchor fixture |
| Final workflow state | Full three-config dry run reports all requested outputs present and up to date |
| Publication check | Exported the staged repository to a fresh temporary workspace, deployed the bundle, generated its fixture and obtained the documented 46-job smoke dry run; Python/JSON/YAML/bash syntax and bundle checksums pass |

The legacy comparison changes one acceptance decision: `759 CA>C` becomes unresolved. All 44 biallelic records have some count changes with stricter classification; other statuses agree. The accepted `560 CCT>C` retains 757 exact HP2 ALT reads, ALT fraction 0.8885 and OTHER fraction 0.1125.

Two separate leave-one-out tagging diagnostics retain the original phase orientation while omitting the assessed site from tagging input:

- `560 CCT>C` remains ACCEPT: withheld HP2 ALT fraction 0.7218; between-HP difference 0.7107.
- `714 G>C`, the weakest accepted SNV by expected ALT fraction, remains ACCEPT: withheld HP2 ALT fraction 0.3841; difference 0.3794.
- All 85,890 alignment identities remain accounted for and production VCF/support hashes stay unchanged. Withholding 560 changes scratch `727 G>A` to ACCEPT, demonstrating assignment sensitivity; its production uncertainty mask is retained. Withholding 714 changes no catalogue status.

## Documentation and validation limits

README updates explain the new support method, phase guard, deletion-mask precedence, test commands and setup bundle. `docs/validation/KEL_indel_revalidation.md` and the bundle runbook/validation records describe the actual code, logs, outputs and reproducible commands. The full feature/change description is retained in `docs/changes/KEL_apple_silicon_features.md`.

These checks establish computational behavior and internal consistency. The synthetic fixtures are small engineered cases, and the leave-one-out diagnostics retain the existing phase solution. External KEL genotype/phase truth, switch-error rate, representative sensitivity/specificity and assay validation remain unavailable. Historical ABO validation is narrative evidence from the repository, not a newly reproduced dataset. Native Clair3 v2/PyTorch behavior is documented without claiming equivalence to earlier TensorFlow-based results.
