import gzip
from collections import defaultdict
from pathlib import Path


def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "rt")


def parse_vcf(path):
    records = []
    with open_text(path) as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt = fields[4]
            alts = [] if alt in ("", ".") else alt.split(",")

            gt = "."
            ps = "."
            if len(fields) >= 10:
                fmt = fields[8].split(":")
                sample = fields[9].split(":")
                values = dict(zip(fmt, sample))
                gt = values.get("GT", ".")
                ps = values.get("PS", ".")

            records.append(
                {
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alts": alts,
                    "gt": gt,
                    "ps": ps,
                }
            )
    return records


def genotype_alleles(gt):
    if not gt or gt == ".":
        return None, None
    sep = "|" if "|" in gt else "/" if "/" in gt else None
    if sep is None:
        return None, None
    parts = gt.split(sep)
    if len(parts) != 2 or any(p == "." for p in parts):
        return None, sep
    try:
        return tuple(int(p) for p in parts), sep
    except ValueError:
        return None, sep


def is_heterozygous(gt):
    alleles, _ = genotype_alleles(gt)
    return alleles is not None and alleles[0] != alleles[1]


def is_homozygous_alt(gt):
    alleles, _ = genotype_alleles(gt)
    return (
        alleles is not None
        and alleles[0] == alleles[1]
        and alleles[0] != 0
    )


def is_phased_heterozygous(gt):
    alleles, sep = genotype_alleles(gt)
    return (
        alleles is not None
        and sep == "|"
        and alleles[0] != alleles[1]
    )


analysis = snakemake.wildcards.analysis
full_records = parse_vcf(snakemake.input.full_vcf)
ready_records = parse_vcf(snakemake.input.phasing_ready_vcf)
phased_records = parse_vcf(snakemake.input.phased_vcf)

full_positions = {(r["chrom"], r["pos"]) for r in full_records}
full_biallelic = sum(len(r["alts"]) == 1 for r in full_records)
full_multiallelic = sum(len(r["alts"]) > 1 for r in full_records)

ready_het = sum(is_heterozygous(r["gt"]) for r in ready_records)
ready_hom_alt = sum(is_homozygous_alt(r["gt"]) for r in ready_records)

phased_het_records = [r for r in phased_records if is_phased_heterozygous(r["gt"])]
unphased_het = sum(
    is_heterozygous(r["gt"]) and not is_phased_heterozygous(r["gt"])
    for r in phased_records
)

phase_sets = defaultdict(list)
for record in phased_het_records:
    if record["ps"] not in ("", "."):
        phase_sets[record["ps"]].append(record)

largest_ps = "."
largest_records = []
if phase_sets:
    largest_ps, largest_records = max(
        phase_sets.items(),
        key=lambda item: (
            len(item[1]),
            max(r["pos"] for r in item[1]) - min(r["pos"] for r in item[1]),
        ),
    )

if largest_records:
    block_start = min(r["pos"] for r in largest_records)
    block_end = max(r["pos"] for r in largest_records)
    block_span = block_end - block_start + 1
else:
    block_start = "."
    block_end = "."
    block_span = "."

phased_fraction = len(phased_het_records) / ready_het if ready_het else 0.0

columns = [
    "analysis",
    "total_variant_records",
    "unique_variant_positions",
    "biallelic_records",
    "multiallelic_records",
    "phasing_ready_records",
    "heterozygous_biallelic_records",
    "homozygous_alt_biallelic_records",
    "phased_heterozygous_records",
    "unphased_heterozygous_records",
    "phased_fraction",
    "number_of_phase_sets",
    "largest_phase_set",
    "largest_phase_block_variants",
    "largest_phase_block_start",
    "largest_phase_block_end",
    "largest_phase_block_span_bp",
]

values = [
    analysis,
    len(full_records),
    len(full_positions),
    full_biallelic,
    full_multiallelic,
    len(ready_records),
    ready_het,
    ready_hom_alt,
    len(phased_het_records),
    unphased_het,
    f"{phased_fraction:.6f}",
    len(phase_sets),
    largest_ps,
    len(largest_records),
    block_start,
    block_end,
    block_span,
]

output = Path(str(snakemake.output.tsv))
output.parent.mkdir(parents=True, exist_ok=True)
with output.open("w") as handle:
    handle.write("\t".join(columns) + "\n")
    handle.write("\t".join(map(str, values)) + "\n")
