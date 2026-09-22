# Synthetic core workflow smoke fixture

This fixture exercises input staging, NanoPlot, minimap2, BAM filtering/indexing,
target generation, alignment/depth QC, and summary generation. It does not validate
biological accuracy, KEL interpretation, Clair3 calling, phasing, or consensus.
All variant calling, phasing, haplotype reconstruction, and consensus flags are
disabled. The core workflow still prepares variant and phasing BAMs for QC.

`generate_fixture.py` uses seed `20260921` to create a 24,000 bp nonrepetitive random
reference and two overlapping 13,000 bp amplicons. Each amplicon contains 24 reads,
with balanced forward/reverse orientations, modest end trimming, deterministic
0.4–0.8% substitutions, and varied per-base qualities. Four full-length reads per
amplicon retain target-edge coverage. Each input directory has two `.fastq` chunks
and one `.fastq.gz` chunk, with 8 reads in every chunk. Gzip timestamps are fixed.

Expected analysis: `SYNTHETIC__TEST`.

| Amplicon | Reference coordinates, 1-based inclusive | Reads | Chunks |
| --- | --- | ---: | ---: |
| amplicon1 | SYNTHETIC_REF:1001-14000 | 24 | 3 |
| amplicon2 | SYNTHETIC_REF:10001-23000 | 24 | 3 |

Expected totals: 48 staged reads, 48 unique primary mapped reads, 48 reads passing
MAPQ 30 filtering, and 48 reads passing the 8,000 bp phasing-length threshold.
The expected reference length is 24,000 bp. The merged target has 22,000 bases:

```text
SYNTHETIC_REF	1000	23000
```

The expected mapping/filtering counts are tool checks to verify after execution;
they are not pre-existing run results. `expected.json` records exact generated
length/base totals; `expected_reads.tsv` records every source span and error count;
`sha256.tsv` records reference and FASTQ hashes. `samples.tsv` uses absolute input
paths, so the workflow may be run in a separate smoke working directory.

Regenerate only this fixture with:

```bash
python3 setup/smoke/generate_fixture.py
```

`config.yaml` sets `qc.nanoplot_no_static: true` for the optional rule flag being
added during setup. NanoPlot statistics and interactive reports remain useful;
static plot export is unnecessary for this execution smoke test.
