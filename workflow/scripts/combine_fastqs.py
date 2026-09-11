#!/usr/bin/env python3
"""Stream one or more FASTQ/FASTQ.GZ inputs into one gzipped FASTQ.

Works both as a Snakemake script and as a standalone CLI, which makes the
input layer easy to unit-test without requiring Snakemake to be installed.
"""

import argparse
import gzip
import shutil
from pathlib import Path


def is_gzip(path):
    return str(path).lower().endswith(".gz")


def combine(inputs, output):
    inputs = [str(x) for x in inputs]
    if not inputs:
        raise ValueError("No FASTQ inputs supplied")

    out = Path(output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(out, "wb", compresslevel=6) as dst:
        for src_name in inputs:
            src = Path(src_name)
            opener = gzip.open if is_gzip(src) else open
            with opener(src, "rb") as handle:
                shutil.copyfileobj(handle, dst, length=1024 * 1024)

    if out.stat().st_size == 0:
        raise RuntimeError(f"Combined FASTQ is empty: {out}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-o", "--output", required=True)
    ap.add_argument("inputs", nargs="+")
    args = ap.parse_args()
    combine(args.inputs, args.output)


if "snakemake" in globals():
    Path(snakemake.log[0]).parent.mkdir(parents=True, exist_ok=True)
    with open(snakemake.log[0], "w") as log:
        log.write(f"Combining {len(snakemake.input.fastqs)} FASTQ chunk(s)\n")
        for p in snakemake.input.fastqs:
            log.write(f"  {p}\n")
    combine(list(snakemake.input.fastqs), snakemake.output.fastq)
elif __name__ == "__main__":
    main()
