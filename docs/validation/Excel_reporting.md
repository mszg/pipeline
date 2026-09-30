# Excel evaluation reporting — implementation verification

**Reporting implementation and artifact checks passed against existing KEL and
ABO evidence.** These checks establish report construction and consistency with
available analysis artifacts. They do not establish biological genotype/phase
truth or validated assay acceptance thresholds.

This document records **prior implementation evidence from 2026-09-29**. Tests
and artifact checks rerun for the versioned release are recorded separately in
[the v0.7 release notes](../releases/v0.7.md). The prior results below must not be
read as independent verification of a later source tree.

## Scope and files

- `Snakefile` and `workflow/rules/report.smk`: strict final `evaluation_report`
  rule, explicit prerequisites, workbook and JSON included in the default target.
- `workflow/scripts/evaluation_metrics.py`: read-only collection from structured
  evidence, intended-target coverage, explicit missing/malformed states and source
  hashes. No BAM/FASTQ reprocessing.
- `workflow/scripts/evaluation_policy.py`: manifest identity, separate sample and
  control policies, gene overrides and transparent status aggregation.
- `workflow/scripts/evaluation_workbook.py`: five styled, filterable sheets,
  literal identifiers/strings, safe links and atomic XLSX writing.
- `workflow/scripts/evaluation_report.py`: Snakemake entry point, explicit
  standalone partial CLI, JSON evidence, and separation of report provenance from
  available analysis provenance.
- `workflow/envs/evaluation_report.yaml` and `config/config.yaml`: reporting
  dependencies and configuration. No validated acceptance cutoff is introduced.
- `test/test_evaluation_{policy,workbook,report}.py` and
  `test/verify_evaluation_workbook.py`: synthetic tests and independent artifact
  checks. The README documents generation commands, partial behavior and policy.

The feature does not change variant calling, filtering, the phase guard,
haplotagging, allele support or consensus behavior. No allele-name or phenotype
interpretation was added.

## Prior implementation checks

The implementation was exercised on Apple Silicon arm64 using an existing
isolated environment with Python 3.12.14 and Snakemake 9.27.0. Reporting libraries
were openpyxl 3.1.5 and et-xmlfile 2.0.0, installed from pure Python wheels; PyYAML
was already available. The dedicated reporting Conda definition declares Python
3.12, openpyxl 3.1.5 and PyYAML >=6.

| Prior check | Recorded result |
| --- | --- |
| Focused reporting tests | 48 passed |
| Entire repository unit-test discovery | 86 passed, including 38 existing analysis tests |
| Strict dry runs on existing KEL and ABO outputs | Each planned one report job; no analysis reruns |
| Strict report execution | Two actual workbooks generated; each completed its report job |
| Final default-target dry runs | Both current, including report outputs |
| Explicit partial CLI | Both produced incomplete reports with execution unverified |
| Independent workbook reopening | All five sheets, expected manifest coverage, row counts, filters, freeze panes, numeric values, literal strings, absence of formulas and links passed |
| Independent metric comparisons | Read counts, coverage, VCF annotations, phasing, assignment, support decisions, haplotype depth, FASTA statistics and reference-coordinate masks matched their sources |
| Metric-source integrity | Recorded source hashes, sizes and modification times matched at verification |
| Preservation against pre-reporting records | All 73 original FASTQs and 58 baseline artifacts/audits retained hashes, sizes and modification times; 23 existing analysis source files retained hashes |

The five checked sheets were `Sample_Summary`, `Gene_QC`, `Variants`,
`Haplotypes` and `Run_Info`. The preservation checks used existing input audits
and an existing source snapshot. `Snakefile` and `workflow/rules/report.smk` were
intentionally changed to integrate the final reporting step; reporting modules,
configuration and documentation were additional changes.

Execution and QC are reported separately. Strict execution was labelled
`DEPENDENCIES SATISFIED (cached outputs may be reused)`, which does not claim that
prior analysis files ran anew. Current KEL/ABO technical QC is **NOT ASSESSED**
because the required minimum `target_breadth_fraction` is unconfigured. The
existing 50× haplotype callability depth does not establish a validated minimum
gene-wide breadth fraction. Unavailable reporting metadata is explicitly labelled
rather than inferred.

Synthetic tests cover zero variants, no heterozygotes, a sample with both ABO
and KEL, a missing expected gene, status precedence with every reason retained,
multiallelic and missing/zero annotations, overlapping targets and zero-depth
positions, truncated depth, duplicate input sources, leading-zero identifiers,
repeated control labels, separate control policies, formula-like external strings,
empty references, missing expected targets, inconsistent haplotype coverage,
assessment attribution, and reference-coordinate masks despite indel-related
FASTA coordinate shifts. Fixture workbooks are labelled `SYNTHETIC DEMONSTRATION`.

## Reproduce checks on local artifacts

From the repository root, in an environment containing the declared reporting
libraries and the dependencies needed by the analysis unit tests:

```bash
python -m unittest discover -s test -p 'test_evaluation*.py' -v
python -m unittest discover -s test -p 'test_*.py' -v

python test/verify_evaluation_workbook.py \
  /path/to/run/results/reports/pipeline_evaluation.xlsx \
  --output /path/to/local/report-verification.json
```

Replace the absolute placeholder paths with the local run and verification
locations. Run the verifier separately for each workbook. It expects a complete
strict report with calling, phasing and reconstruction enabled; partial reports
and disabled-stage behavior are covered by synthetic tests. Generation and
regeneration commands are in the README.

Real workbooks, JSON evidence, raw sequencing files and local verification logs
are not release assets. This public record intentionally omits individual
sample results and genotype examples.

## Limits

Workbooks were reopened programmatically with openpyxl and their contents and
formatting metadata checked. No interactive desktop Excel or LibreOffice
inspection was performed. No fresh solve/build of the dedicated per-rule Conda
environment was performed; execution used the existing isolated environment with
the reporting libraries installed. Real multi-sample/control reporting was not
tested; those cases used synthetic fixtures.

Unique-molecule counts, haplotype median depth and original analysis commit
identity are unavailable in existing standard structured artifacts. Only tool
provenance actually recorded in VCF headers is presented as analysis provenance;
report-environment versions are explicitly separate. Missing annotations are not
fabricated, and available depth does not establish independent molecule support.

Reconstruction still requires an informative phase block under the existing
pipeline guard: exactly one informative phase set per target contig. The report
neither relaxes that guard nor treats absence of heterozygotes as an automatic QC
failure. Explicit partial reporting remains available when upstream prerequisites
are unavailable. No independent genotype/phase truth, clinical validation or
validated assay acceptance thresholds are established by these reporting checks.
