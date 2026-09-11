import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_utils import parse_flagstat, parse_coverage

amp_rows = []
with open(snakemake.input.amplicon_summary, newline="") as handle:
    amp_rows = list(csv.DictReader(handle, delimiter="\t"))

analyses = snakemake.params.analyses.split(",") if snakemake.params.analyses else []

outpath = Path(snakemake.output.tsv)
outpath.parent.mkdir(parents=True, exist_ok=True)
with outpath.open("w", newline="") as out:
    fields = [
        "sample", "gene", "analysis", "amplicons", "total_input_reads",
        "mapped_percent", "reference_coverage_percent", "mean_depth"
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for i, analysis in enumerate(analyses):
        subset = [r for r in amp_rows if f'{r["sample"]}__{r["gene"]}' == analysis]
        if subset:
            sample, gene = subset[0]["sample"], subset[0]["gene"]
        else:
            sample, gene = analysis, "NA"
        numeric_reads = []
        for r in subset:
            try:
                numeric_reads.append(int(float(r["reads"])))
            except (ValueError, TypeError):
                pass
        _, _, mapped_pct = parse_flagstat(snakemake.input.flagstats[i])
        coverage_pct, mean_depth = parse_coverage(snakemake.input.coverages[i])
        writer.writerow({
            "sample": sample,
            "gene": gene,
            "analysis": analysis,
            "amplicons": len(subset),
            "total_input_reads": sum(numeric_reads) if numeric_reads else "NA",
            "mapped_percent": mapped_pct,
            "reference_coverage_percent": coverage_pct,
            "mean_depth": mean_depth,
        })
