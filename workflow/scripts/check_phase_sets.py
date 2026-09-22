"""Require one connected informative phase set per target contig before HP splitting."""
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def check_phase_sets(vcf_path, target_bed, expected_sample=None):
    targets = {line.split()[0] for line in Path(target_bed).read_text().splitlines() if line.strip() and not line.startswith("#")}
    if not targets:
        raise ValueError("Haplotype reconstruction requires nonempty target intervals")
    phases = defaultdict(lambda: defaultdict(list))
    unphased = defaultdict(int)
    samples = None
    opener = gzip.open if str(vcf_path).endswith(".gz") else open
    with opener(vcf_path, "rt") as handle:
        for line in handle:
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if len(samples) != 1 or (expected_sample and samples[0] != expected_sample):
                    raise ValueError(f"Expected one VCF sample {expected_sample or ''!r}; found {samples}")
            elif not line.startswith("#") and line.strip():
                fields = line.rstrip().split("\t")
                if samples is None or len(fields) != 10:
                    raise ValueError("Malformed single-sample phased VCF")
                chrom, position = fields[0], int(fields[1])
                if chrom not in targets:
                    raise ValueError(f"VCF contig {chrom} is absent from the target BED")
                values = dict(zip(fields[8].split(":"), fields[9].split(":")))
                gt = values.get("GT", ".")
                alleles = gt.replace("|", "/").split("/")
                if len(alleles) != 2 or "." in alleles or alleles[0] == alleles[1]:
                    continue
                if "|" not in gt:
                    unphased[chrom] += 1
                    continue
                ps = values.get("PS", ".")
                if ps in ("", "."):
                    raise ValueError(f"Phased heterozygote {chrom}:{position} has no PS")
                phases[chrom][ps].append(position)
    if samples is None:
        raise ValueError("VCF sample header is missing")
    rows = []
    for chrom in sorted(targets):
        groups = phases[chrom]
        if len(groups) != 1:
            raise ValueError(
                f"Cannot reconstruct consistent HP1/HP2 for {chrom}: expected one informative "
                f"phase set, found {len(groups)} ({', '.join(sorted(groups)) or 'none'}). "
                "Disconnected blocks require separate reconstruction."
            )
        ps, positions = next(iter(groups.items()))
        rows.append({"contig": chrom, "phase_set": ps, "phased_heterozygotes": len(positions),
                     "unphased_heterozygotes": unphased[chrom], "start": min(positions),
                     "end": max(positions), "status": "PASS"})
    return rows


if "snakemake" in globals():
    rows = check_phase_sets(snakemake.input.vcf, snakemake.input.bed, snakemake.wildcards.analysis)
    output = Path(str(snakemake.output.tsv))
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
