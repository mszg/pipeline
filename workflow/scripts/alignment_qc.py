#!/usr/bin/env python3
"""Alignment diagnostics tailored to long-range ONT amplicons."""

import csv
from pathlib import Path
import pysam

bam = str(snakemake.input.bam)
summary_out = Path(str(snakemake.output.summary))
length_out = Path(str(snakemake.output.length_bins))
summary_out.parent.mkdir(parents=True, exist_ok=True)
length_out.parent.mkdir(parents=True, exist_ok=True)

bins = [
    ("<2kb", 0, 2000),
    ("2-5kb", 2000, 5000),
    ("5-8kb", 5000, 8000),
    ("8-11kb", 8000, 11000),
    ("11-13kb", 11000, 13000),
    (">=13kb", 13000, None),
]
length_counts = {name: [0, 0] for name, _, _ in bins}  # primary mapped, MAPQ>=30

stats = {
    "total_alignment_records": 0,
    "primary_records": 0,
    "primary_mapped": 0,
    "mapq_ge_20": 0,
    "mapq_ge_30": 0,
    "mapq_ge_50": 0,
    "secondary": 0,
    "supplementary": 0,
}

high_lengths = []
low_lengths = []

with pysam.AlignmentFile(bam, "rb") as handle:
    for read in handle.fetch(until_eof=True):
        stats["total_alignment_records"] += 1
        if read.is_secondary:
            stats["secondary"] += 1
        if read.is_supplementary:
            stats["supplementary"] += 1

        if read.is_secondary or read.is_supplementary:
            continue
        stats["primary_records"] += 1
        if read.is_unmapped:
            continue
        stats["primary_mapped"] += 1

        mq = int(read.mapping_quality)
        if mq >= 20:
            stats["mapq_ge_20"] += 1
        if mq >= 30:
            stats["mapq_ge_30"] += 1
        if mq >= 50:
            stats["mapq_ge_50"] += 1

        qlen = int(read.query_length or 0)
        if mq >= 30:
            high_lengths.append(qlen)
        if mq < 20:
            low_lengths.append(qlen)

        for name, low, high in bins:
            if qlen >= low and (high is None or qlen < high):
                length_counts[name][0] += 1
                if mq >= 30:
                    length_counts[name][1] += 1
                break

stats["mean_length_mapq_ge_30"] = (
    f"{sum(high_lengths) / len(high_lengths):.2f}" if high_lengths else "NA"
)
stats["mean_length_mapq_lt_20"] = (
    f"{sum(low_lengths) / len(low_lengths):.2f}" if low_lengths else "NA"
)

with summary_out.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["metric", "value"])
    for key, value in stats.items():
        writer.writerow([key, value])

with length_out.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["length_bin", "primary_mapped_reads", "mapq_ge_30", "mapq_ge_30_percent"])
    for name, _, _ in bins:
        total, high = length_counts[name]
        pct = 100.0 * high / total if total else 0.0
        writer.writerow([name, total, high, f"{pct:.2f}"])
