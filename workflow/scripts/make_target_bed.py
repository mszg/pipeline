#!/usr/bin/env python3
"""Create a merged BED calling target plus an auditable per-amplicon manifest.

Input target regions are stored in samples.tsv as 1-based inclusive
CONTIG:START-END coordinates. Clair3 expects BED regions, therefore this script
converts them to 0-based half-open intervals and merges overlapping/adjacent
intervals on the same contig.
"""

import csv
import re
from pathlib import Path

REGION_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")

samples_tsv = Path(str(snakemake.input.samples))
analysis = str(snakemake.params.analysis)
bed_path = Path(str(snakemake.output.bed))
manifest_path = Path(str(snakemake.output.manifest))

records = []
with samples_tsv.open(newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    for row in reader:
        row_analysis = f"{row.get('sample', '').strip()}__{row.get('gene', '').strip()}"
        if row_analysis != analysis:
            continue
        region = row.get("target_region", "").strip()
        if not region:
            continue
        match = REGION_RE.match(region)
        if not match:
            raise ValueError(
                f"Invalid target_region={region!r} for {analysis}; expected CONTIG:START-END"
            )
        contig, start_s, end_s = match.groups()
        start, end = int(start_s), int(end_s)
        if start < 1 or end < start:
            raise ValueError(f"Invalid target_region coordinates: {region}")
        records.append(
            {
                "sample": row["sample"].strip(),
                "gene": row["gene"].strip(),
                "amplicon": row["amplicon"].strip(),
                "region": region,
                "contig": contig,
                "start_1based": start,
                "end_1based": end,
                "start0": start - 1,
                "end0": end,
            }
        )

bed_path.parent.mkdir(parents=True, exist_ok=True)
manifest_path.parent.mkdir(parents=True, exist_ok=True)

with manifest_path.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t", lineterminator="\n")
    writer.writerow(
        [
            "sample",
            "gene",
            "amplicon",
            "target_region",
            "contig",
            "start_1based_inclusive",
            "end_1based_inclusive",
            "bed_start_0based",
            "bed_end_half_open",
        ]
    )
    for r in records:
        writer.writerow(
            [
                r["sample"],
                r["gene"],
                r["amplicon"],
                r["region"],
                r["contig"],
                r["start_1based"],
                r["end_1based"],
                r["start0"],
                r["end0"],
            ]
        )

# Merge overlapping or directly adjacent BED intervals by contig.
merged = []
for contig in sorted({r["contig"] for r in records}):
    intervals = sorted(
        (r["start0"], r["end0"]) for r in records if r["contig"] == contig
    )
    for start0, end0 in intervals:
        if not merged or merged[-1][0] != contig or start0 > merged[-1][2]:
            merged.append([contig, start0, end0])
        else:
            merged[-1][2] = max(merged[-1][2], end0)

with bed_path.open("w") as out:
    for contig, start0, end0 in merged:
        out.write(f"{contig}\t{start0}\t{end0}\n")
