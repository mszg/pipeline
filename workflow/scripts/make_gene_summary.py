import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_utils import parse_flagstat, parse_coverage

with open(snakemake.input.amplicon_summary, newline="") as handle:
    amp_rows = list(csv.DictReader(handle, delimiter="\t"))

analyses = snakemake.params.analyses.split(",") if snakemake.params.analyses else []


def as_int(value):
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return None


def pct(num, den):
    if num is None or den in (None, 0):
        return "NA"
    return f"{100.0 * num / den:.2f}"


outpath = Path(snakemake.output.tsv)
outpath.parent.mkdir(parents=True, exist_ok=True)
with outpath.open("w", newline="") as out:
    fields = [
        "sample", "gene", "analysis", "amplicons", "total_input_reads",
        "primary_mapped_reads", "primary_mapped_percent",
        "variant_analysis_reads", "phasing_reads",
        "raw_alignment_mapped_percent",
        "raw_reference_coverage_percent", "raw_mean_depth",
        "analysis_reference_coverage_percent", "analysis_mean_depth",
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for i, analysis in enumerate(analyses):
        subset = [r for r in amp_rows if f'{r["sample"]}__{r["gene"]}' == analysis]
        if subset:
            sample, gene = subset[0]["sample"], subset[0]["gene"]
        else:
            sample, gene = analysis, "NA"

        input_reads = [as_int(r.get("reads")) for r in subset]
        input_reads = [x for x in input_reads if x is not None]
        primary_reads = [as_int(r.get("primary_mapped_reads")) for r in subset]
        primary_reads = [x for x in primary_reads if x is not None]
        total_input = sum(input_reads) if input_reads else None
        total_primary = sum(primary_reads) if primary_reads else None

        _, _, raw_mapped_pct = parse_flagstat(snakemake.input.raw_flagstats[i])
        raw_cov_pct, raw_depth = parse_coverage(snakemake.input.raw_coverages[i])
        variant_total, _, _ = parse_flagstat(snakemake.input.variant_flagstats[i])
        phasing_total, _, _ = parse_flagstat(snakemake.input.phasing_flagstats[i])
        analysis_cov_pct, analysis_depth = parse_coverage(snakemake.input.variant_coverages[i])

        writer.writerow({
            "sample": sample,
            "gene": gene,
            "analysis": analysis,
            "amplicons": len(subset),
            "total_input_reads": total_input if total_input is not None else "NA",
            "primary_mapped_reads": total_primary if total_primary is not None else "NA",
            "primary_mapped_percent": pct(total_primary, total_input),
            "variant_analysis_reads": variant_total,
            "phasing_reads": phasing_total,
            "raw_alignment_mapped_percent": raw_mapped_pct,
            "raw_reference_coverage_percent": raw_cov_pct,
            "raw_mean_depth": raw_depth,
            "analysis_reference_coverage_percent": analysis_cov_pct,
            "analysis_mean_depth": analysis_depth,
        })
