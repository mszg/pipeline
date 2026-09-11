import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_utils import parse_nanostats, parse_flagstat, parse_coverage

rows = []
with open(snakemake.params.samples_tsv, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("sample", "").strip():
            row["analysis"] = f'{row["sample"]}__{row["gene"]}'
            row["unit"] = f'{row["analysis"]}__{row["amplicon"]}'
            rows.append(row)

rows_by_unit = {r["unit"]: r for r in rows}
units = snakemake.params.units.split(",") if snakemake.params.units else []

chunk_counts = {}
with open(snakemake.input.manifest, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        chunk_counts[row["unit"]] = chunk_counts.get(row["unit"], 0) + 1

outpath = Path(snakemake.output.tsv)
outpath.parent.mkdir(parents=True, exist_ok=True)
with outpath.open("w", newline="") as out:
    fields = [
        "sample", "gene", "amplicon", "unit", "fastq_chunks", "reads",
        "mean_read_length", "n50", "mean_read_quality", "mapped_percent",
        "reference_coverage_percent", "mean_depth"
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for i, unit in enumerate(units):
        meta = rows_by_unit[unit]
        ns = parse_nanostats(snakemake.input.nanostats[i])
        _, _, mapped_pct = parse_flagstat(snakemake.input.flagstats[i])
        coverage_pct, mean_depth = parse_coverage(snakemake.input.coverages[i])
        writer.writerow({
            "sample": meta["sample"],
            "gene": meta["gene"],
            "amplicon": meta["amplicon"],
            "unit": unit,
            "fastq_chunks": chunk_counts.get(unit, 0),
            "reads": ns["reads"],
            "mean_read_length": ns["mean_read_length"],
            "n50": ns["n50"],
            "mean_read_quality": ns["mean_read_quality"],
            "mapped_percent": mapped_pct,
            "reference_coverage_percent": coverage_pct,
            "mean_depth": mean_depth,
        })
