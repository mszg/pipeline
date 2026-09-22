# KEL indel repair and reconstruction validation — 2026-09-21

This record describes observed local results after repair of the v0.6 counter.
It establishes regression behavior, workflow execution and internal consistency;
there is no orthogonal KEL truth set for variant or switch-error accuracy.
Historical ABO outputs were not available in this clone and were not reproduced.

## Inputs and configuration

One biological sample, with technical analysis key `KEL__KEL`; no patient/sample
identifier was assigned. All 36 original FASTQs (399,103 reads) remain unchanged
against their pre-run SHA-256 inventory. The reference is NCBI RefSeqGene
`NG_007492.3`, 28,313 bp. Conservative outer targets are `411-15302` and
`14079-28023`; merged BED is `NG_007492.3 410 28023`. Coordinates are RefSeqGene.
The ONT model is `dna_r10.4.1_e8.2_400bps_hac@v5.2.0`.

Native Apple Silicon core tools and Clair3 2.0.3 run in separate isolated Conda
environments. The local configs are applied in this order: `kel.yaml`,
`kel.calling.yaml`, `kel.reconstruction.yaml`. Existing thresholds are unchanged:
50x callable depth, 20 exact REF+ALT observations per HP, 0.30 expected heterozygous
ALT fraction, 0.25 between-HP ALT-fraction difference, 0.80 homozygous ALT fraction,
at most 0.25 OTHER fraction, and caller PASS required.

## Repairs and tests

- Replaced anchor-only indel counting with exact aligned allele windows spanning
  anchor, indel/repeat context and right flank. Incomplete windows, wrong anchors,
  ambiguous bases and competing sequence are OTHER. Equivalent shifted indels in
  repeats cannot be mistaken for exact REF. Only normalized pure small indels are
  supported; neighboring variants can conservatively reduce exact support.
- Tightened SNV pileup handling so placeholders/skips, ambiguous bases and
  competing indel events cannot be exact REF/ALT.
- Added a fail-fast check for exactly one informative phase set per target contig
  before HP splitting. HP labels do not establish phase between different contigs.
- Corrected masking of accepted deletions with zero aligned base depth. Rescue
  requires enough exact ALT reads on that haplotype, a callable anchor and no
  immutable mask over the REF span. Blocked deletions retain their low-depth mask.

All **38 unit tests passed**, covering read boundaries, skips, placeholders,
repeat-shifted indels, long events, reference mismatch, primary read eligibility,
phase blocks and mask precedence. The portable integration test runs the actual
support, coverage and mask scripts and bcftools consensus against known synthetic
alleles. Both full 500 bp haplotype sequences matched expected truth. A separate
blocked-anchor deletion fixture also matched the fully masked expected allele.
The synthetic test supplies known HP tags; actual haplotag execution is covered
by the KEL run below. These small fixtures are not a representative ONT benchmark.

```bash
python -m unittest discover -s test -p 'test_*.py' -v
python test/integration_haplotype_consensus.py --output /path/to/empty/directory
```

## Real KEL results

The reconstruction dry run planned 11 jobs. All completed successfully in
53.8 seconds; after tightening mask precedence, the five affected final jobs
were dry-run and completed again. Calling/phasing inputs were retained.

| Metric | Observed result |
| --- | --- |
| Normalized catalogue | 45 records: 44 biallelic, 1 multiallelic |
| Phasing | 41/41 eligible heterozygotes in PS516, positions 516–25592 |
| Input phasing reads | 85,890 |
| HP1 / HP2 / unassigned reads | 32,100 / 21,523 / 32,267 |
| Support decisions | 33 ACCEPT: 32 SNVs and one deletion; 12 UNRESOLVED |
| HP1 minimum / mean target depth | 630 / 13,371.93 |
| HP2 minimum / mean target depth | 13 / 8,762.21 |
| Target positions below 50x | HP1: 0; HP2: 87 (`411-497`) |
| HP1 FASTA | 28,313 bases; 724 Ns, of which 24 target positions |
| HP2 FASTA | 28,311 bases; 811 Ns, of which 111 target positions |

Each FASTA includes 700 outside-target masked reference positions. HP2 includes
the accepted two-base deletion `560 CCT>C`; all accepted ALT variants were applied
on their designated haplotypes. No KEL deletion required the zero-depth rescue;
that defect is exercised by the synthetic fixture. Unresolved positions are
624, 639, 660, 692, 696, 717, 727, 758, 759, 1688, 25592 and 27258. The last is
multiallelic; the other eleven fail current read-support criteria. The uncertainty
mask covers 24 reference bases, and is retained in both FASTAs. Inserted sequence
ambiguity remains in the catalogue/support table rather than being encoded by a
reference-coordinate mask.

Re-running the archived defective counter on identical HP BAMs changed one
decision: `759 CA>C`, formerly ACCEPT, is now UNRESOLVED. HP2 REF/ALT/OTHER changed
from `99/860/1` to `4/4/952`. The full allele window is inconsistent in nearly all
HP2 reads, so those observations no longer support an exact deletion allele.
This is a conservative support decision, not proof that the underlying biological
deletion is absent. The accepted `560 CCT>C` retains HP2 exact ALT support of 757
reads, ALT fraction 0.8885, and OTHER fraction 0.1125. No acceptance threshold was
relaxed. All 44 biallelic records have some count changes; other statuses agree.

Independent output checks passed: read partitions account for all 85,890 input
reads; assigned reads retain PS516; HP BAM indexes and VCF indexes parse; all
catalogue records are accounted for; accepted VCF alleles/genotypes match the
support table; uncertainty and final masks match their intended sets; and both
complete FASTA sequences equal an independent reference-coordinate reconstruction.
Indel denominators were checked for every site/HP, and 200 read windows were also
checked with an independent aligned-pair extraction. Original FASTQ hashes were
verified again after reconstruction.

Two additional leave-one-out diagnostics separately removed `560 CCT>C` and the
weakest accepted SNV by expected ALT fraction, `714 G>C`, from the haplotag input.
Evaluation retained the original phased genotype. Both remain ACCEPT: withheld
HP2 ALT fractions are 0.7218 and 0.3841, with HP differences 0.7107 and 0.3794.
All 85,890 alignment identities and PS516 were preserved. Withholding 560 changed
scratch `727 G>A` from UNRESOLVED to ACCEPT, showing assignment sensitivity at
that site; its production mask remains unchanged. Withholding 714 changed no
catalogue status. These diagnostics are conditional on the original phase
solution and do not independently validate phase. Production VCF/support hashes
were unchanged. Detailed evidence is in `setup/indel_fix/leave_one_out/`.

## Evidence and remaining scope

The [published setup bundle](../../deployment/kel-apple-silicon/README.md)
contains the reusable scripts, exact native environment locks and selected
logs/JSON evidence below. Its manifest records publication edits to machine
paths; sequencing files, installed environments and full input audits stay local.

Local evidence is stored outside the clone in workspace `setup/indel_fix/`:
`unit-tests.log`, `portable_synthetic/verification.json`,
`kel-reconstruction-dry-run.log`, `kel-reconstruction-run.log`,
`kel-mask-final-run.log`, `kel_legacy_comparison.json`, and
`kel_reconstruction_verification.json`. The final `kel-final-dry-run.log`
reports nothing to do. Archived original scripts and the
pre-fix failing synthetic result are retained there. Workspace `setup/README.md`
contains commands, exact environment locks and reference/input provenance.

The FASTAs are reference-guided reconstructions with explicit uncertainty, not
finished de novo assemblies. Phase/support uses the same sequencing experiment;
one phase block does not measure switch-error rate. Independent genotype/phase
truth, representative sensitivity/specificity and assay validation remain outside
the evidence available here.
