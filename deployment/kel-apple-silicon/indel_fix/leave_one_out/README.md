# KEL leave-one-out haplotag diagnostics

Run from the workspace root:

```sh
setup/envs/kel-native/bin/python setup/indel_fix/leave_one_out.py > setup/indel_fix/leave_one_out/run.log 2>&1
```

Each case removes one exact variant from a scratch copy of the phased VCF, keeps every other record and its phase orientation, and haplotags the original 85,890-read phasing BAM again. The corrected support counter then evaluates all variants using the original full and phased VCFs and the scratch HP BAMs. Production files and original sequencing files are not edited. Production VCF/support SHA256 checks remained unchanged.

Both withheld variants remained `ACCEPT` under the unchanged production thresholds:

| Withheld variant | HP2 ALT fraction, original → withheld | HP1 ALT fraction, original → withheld | ALT-fraction delta, withheld | HP2 exact ALT / REF / OTHER, withheld |
| --- | --- | --- | --- | --- |
| NG_007492.3:560 CCT>C | 0.8885 → 0.7218 | 0.0036 → 0.0111 | 0.7107 | 755 / 291 / 155 |
| NG_007492.3:714 G>C | 0.5091 → 0.3841 | 0.0044 → 0.0047 | 0.3794 | 474 / 760 / 30 |

All 85,890 alignments remained accounted for. The original HP1 / HP2 / unassigned counts were 32,100 / 21,523 / 32,267. After withholding 560 they were 32,086 / 21,764 / 32,040; after withholding 714 they were 32,029 / 21,826 / 32,035. No read switched directly between HP1 and HP2; assignment changes involved unassigned reads. All tagged reads retained phase set 516.

Withholding 560 also changed the scratch classification of 727 G>A from `UNRESOLVED` to `ACCEPT`: HP2 ALT fraction rose from 0.1566 to 0.3186 and the fraction difference from 0.0951 to 0.2631. This demonstrates sensitivity to the tagging input at 727. Its production classification and mask remain unchanged. Withholding 714 caused no classification changes among any of the 45 catalog variants.

These diagnostics test dependence on direct use of a site during tagging. They remain conditional on the original genotype and phase solution; they do not establish independent phase truth, clinical accuracy, or genotype concordance with an external truth set.

`verification.json` contains both case summaries, assignment transition counts, unchanged-file hashes, and thresholds. Each position directory contains its scratch VCF, tagged and split BAMs, full support tables, exact withheld VCF record, and `commands.log`. `run.log` records the completed execution.
