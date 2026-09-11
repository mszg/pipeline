#!/usr/bin/env python3
"""Create a clean analysis BAM from a coordinate-sorted alignment BAM.

Filtering is intentionally BAM-level so the raw FASTQ and raw BAM remain
untouched and every downstream filtering choice is reproducible.
"""

from pathlib import Path
import re
import pysam

REGION_RE = re.compile(r"^([^:\s]+):(\d+)-(\d+)$")


def parse_region(text):
    text = (text or "").strip()
    if not text:
        return None
    match = REGION_RE.match(text)
    if not match:
        raise ValueError(
            f"Invalid target_region={text!r}; expected CONTIG:START-END (1-based, inclusive)"
        )
    contig, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    if start < 1 or end < start:
        raise ValueError(f"Invalid target_region coordinates: {text}")
    return contig, start - 1, end  # convert to 0-based half-open interval


def overlaps_target(read, region):
    if region is None:
        return True
    contig, start0, end0 = region
    if read.reference_name != contig or read.reference_end is None:
        return False
    return read.reference_start < end0 and read.reference_end > start0


input_bam = str(snakemake.input.bam)
output_bam = str(snakemake.output.bam)
min_mapq = int(snakemake.params.min_mapq)
min_read_length = int(snakemake.params.min_read_length)
target_region = parse_region(str(snakemake.params.target_region))
allow_empty = bool(snakemake.params.allow_empty)

Path(output_bam).parent.mkdir(parents=True, exist_ok=True)
Path(str(snakemake.log[0])).parent.mkdir(parents=True, exist_ok=True)

kept = 0
seen = 0
with pysam.AlignmentFile(input_bam, "rb") as src, pysam.AlignmentFile(
    output_bam, "wb", header=src.header
) as dst:
    for read in src.fetch(until_eof=True):
        seen += 1
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        qlen = read.query_length or 0
        if qlen < min_read_length:
            continue
        if not overlaps_target(read, target_region):
            continue
        dst.write(read)
        kept += 1

if kept == 0 and not allow_empty:
    raise RuntimeError(
        f"Filtering removed every alignment from {input_bam}. "
        f"min_mapq={min_mapq}, min_read_length={min_read_length}, "
        f"target_region={snakemake.params.target_region!r}"
    )

with open(str(snakemake.log[0]), "w") as log:
    log.write(f"input_records\t{seen}\n")
    log.write(f"kept_records\t{kept}\n")
    log.write(f"min_mapq\t{min_mapq}\n")
    log.write(f"min_read_length\t{min_read_length}\n")
    log.write(f"target_region\t{snakemake.params.target_region or 'NONE'}\n")
    log.write(f"allow_empty\t{allow_empty}\n")
