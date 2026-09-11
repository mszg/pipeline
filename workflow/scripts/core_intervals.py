#!/usr/bin/env python3
"""Convert samtools depth output into contiguous high-depth intervals."""

import csv
from pathlib import Path

threshold = int(snakemake.params.threshold)
source = Path(str(snakemake.input.depth))
outpath = Path(str(snakemake.output.tsv))
outpath.parent.mkdir(parents=True, exist_ok=True)

intervals = []
current_contig = None
start = None
end = None


def close_interval():
    global start, end, current_contig
    if start is not None:
        intervals.append((current_contig, start, end, end - start + 1, threshold))
    start = None
    end = None

for line in source.read_text(errors="replace").splitlines():
    if not line:
        continue
    contig, pos_s, depth_s = line.split("\t")[:3]
    pos, depth = int(pos_s), int(depth_s)
    if current_contig is not None and contig != current_contig:
        close_interval()
    current_contig = contig
    if depth >= threshold:
        if start is None:
            start = pos
        end = pos
    else:
        close_interval()
close_interval()

with outpath.open("w", newline="") as out:
    writer = csv.writer(out, delimiter="\t")
    writer.writerow(["contig", "start", "end", "length", "depth_threshold"])
    writer.writerows(intervals)
