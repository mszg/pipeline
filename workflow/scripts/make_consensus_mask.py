from collections import defaultdict
import csv
from pathlib import Path

fai = Path(str(snakemake.input.fai))
target_bed = Path(str(snakemake.input.target_bed))
uncertain_bed = Path(str(snakemake.input.uncertain_bed))
low_depth_bed = Path(str(snakemake.input.low_depth_bed))
support_tsv = Path(str(snakemake.input.support_tsv))
haplotype = int(snakemake.params.haplotype)
callable_min_depth = int(snakemake.params.callable_min_depth)
out = Path(str(snakemake.output.bed))
out.parent.mkdir(parents=True, exist_ok=True)


def read_intervals(path):
    intervals = []
    if not path.exists():
        return intervals
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            intervals.append((fields[0], int(fields[1]), int(fields[2])))
    return intervals


def merge(intervals):
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        if end > start:
            by_chrom[chrom].append((start, end))
    out_intervals = []
    for chrom, values in by_chrom.items():
        current = None
        for start, end in sorted(values):
            if current is None or start > current[1]:
                if current is not None:
                    out_intervals.append((chrom, current[0], current[1]))
                current = [start, end]
            else:
                current[1] = max(current[1], end)
        if current is not None:
            out_intervals.append((chrom, current[0], current[1]))
    return out_intervals


def subtract(intervals, exclusions):
    """Subtract supported deletion spans from low-base-depth intervals only."""
    result = []
    for chrom, start, end in intervals:
        parts = [(start, end)]
        for other_chrom, left, right in exclusions:
            if other_chrom != chrom:
                continue
            updated = []
            for a, b in parts:
                if right <= a or left >= b:
                    updated.append((a, b))
                else:
                    if a < left:
                        updated.append((a, left))
                    if right < b:
                        updated.append((right, b))
            parts = updated
        result.extend((chrom, a, b) for a, b in parts)
    return result


supported_deletions = []
with support_tsv.open() as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row["STATUS"] != "ACCEPT" or row["TYPE"] != "DEL":
            continue
        alleles = row["PHASED_GT"].replace("|", "/").split("/")
        if len(alleles) != 2 or alleles[haplotype - 1] != "1":
            continue
        if int(row[f"HP{haplotype}_ALT"]) < callable_min_depth:
            continue
        ref_allele, alt_allele = row["REF"], row["ALT"]
        if len(alt_allele) != 1 or not ref_allele.startswith(alt_allele):
            raise ValueError("Accepted deletion is not a normalized anchored deletion")
        start = int(row["POS"]) - 1
        # Keep the anchor as well: bcftools skips the entire variant if any
        # reference base in its allele span intersects the consensus mask.
        supported_deletions.append((row["CHROM"], start, start + len(ref_allele)))


contigs = []
lengths = {}
with fai.open() as handle:
    for line in handle:
        fields = line.rstrip("\n").split("\t")
        chrom = fields[0]
        length = int(fields[1])
        contigs.append(chrom)
        lengths[chrom] = length

targets = merge(read_intervals(target_bed))
targets_by_chrom = defaultdict(list)
for chrom, start, end in targets:
    if chrom not in lengths:
        raise ValueError(f"Target BED contig {chrom!r} is not present in reference .fai")
    targets_by_chrom[chrom].append((max(0, start), min(lengths[chrom], end)))

mask = []
for chrom in contigs:
    length = lengths[chrom]
    cursor = 0
    for start, end in targets_by_chrom.get(chrom, []):
        if start > cursor:
            mask.append((chrom, cursor, start))
        cursor = max(cursor, end)
    if cursor < length:
        mask.append((chrom, cursor, length))

mask.extend(read_intervals(uncertain_bed))
# Deletions have no aligned query base at the deleted reference positions.
# Rescue only an accepted ALT on this haplotype with enough exact ALT reads.
# If an immutable mask or low-depth anchor blocks application of the deletion,
# retain its low-depth mask so skipped variants cannot expose reference bases.
# Outside-target and uncertainty masks are never relaxed.
low_depth = read_intervals(low_depth_bed)


def overlaps(chrom, start, end, intervals):
    return any(chrom == other_chrom and start < right and end > left
               for other_chrom, left, right in intervals)


rescuable_deletions = []
for chrom, start, end in supported_deletions:
    if overlaps(chrom, start, end, mask) or overlaps(chrom, start, start + 1, low_depth):
        continue
    rescuable_deletions.append((chrom, start + 1, end))
mask.extend(subtract(low_depth, rescuable_deletions))
mask = merge(mask)

order = {chrom: i for i, chrom in enumerate(contigs)}
mask.sort(key=lambda x: (order.get(x[0], len(order)), x[1], x[2]))

with out.open("w") as handle:
    for chrom, start, end in mask:
        handle.write(f"{chrom}\t{start}\t{end}\n")
