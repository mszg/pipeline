import csv
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from summary_utils import (
    parse_nanostats,
    parse_flagstat,
    parse_coverage,
    parse_metric_tsv,
    longest_interval,
)

rows = []
with open(snakemake.params.samples_tsv, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        if row.get("sample", "").strip():
            row = {k: (v.strip() if isinstance(v, str) else v) for k, v in row.items()}
            row["analysis"] = f'{row["sample"]}__{row["gene"]}'
            row["unit"] = f'{row["analysis"]}__{row["amplicon"]}'
            rows.append(row)

rows_by_unit = {r["unit"]: r for r in rows}
units = snakemake.params.units.split(",") if snakemake.params.units else []

chunk_counts = {}
with open(snakemake.input.manifest, newline="") as handle:
    for row in csv.DictReader(handle, delimiter="\t"):
        chunk_counts[row["unit"]] = chunk_counts.get(row["unit"], 0) + 1


def pct(num, den):
    try:
        num = float(num)
        den = float(den)
        return f"{100.0 * num / den:.2f}" if den else "NA"
    except (TypeError, ValueError):
        return "NA"


outpath = Path(snakemake.output.tsv)
outpath.parent.mkdir(parents=True, exist_ok=True)
with outpath.open("w", newline="") as out:
    fields = [
        "sample", "gene", "amplicon", "unit", "fastq_chunks", "reads",
        "mean_read_length", "n50", "mean_read_quality",
        "primary_mapped_reads", "primary_mapped_percent",
        "mapq_ge_20_reads", "mapq_ge_20_of_primary_percent",
        "mapq_ge_30_reads", "mapq_ge_30_of_primary_percent",
        "mapq_ge_50_reads", "mapq_ge_50_of_primary_percent",
        "secondary_alignments", "supplementary_alignments",
        "mean_length_mapq_ge_30", "mean_length_mapq_lt_20",
        "raw_alignment_mapped_percent",
        "raw_reference_coverage_percent", "raw_mean_depth",
        "analysis_reference_coverage_percent", "analysis_mean_depth",
        "observed_core_interval", "observed_core_length", "core_depth_threshold",
        "target_region", "phasing_min_length",
    ]
    writer = csv.DictWriter(out, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    for i, unit in enumerate(units):
        meta = rows_by_unit[unit]
        ns = parse_nanostats(snakemake.input.nanostats[i])
        _, _, raw_mapped_pct = parse_flagstat(snakemake.input.flagstats[i])
        raw_cov_pct, raw_depth = parse_coverage(snakemake.input.raw_coverages[i])
        analysis_cov_pct, analysis_depth = parse_coverage(snakemake.input.analysis_coverages[i])
        aqc = parse_metric_tsv(snakemake.input.alignment_qc[i])
        core_region, core_length, core_threshold = longest_interval(snakemake.input.core_intervals[i])

        primary = aqc.get("primary_mapped", "NA")
        q20 = aqc.get("mapq_ge_20", "NA")
        q30 = aqc.get("mapq_ge_30", "NA")
        q50 = aqc.get("mapq_ge_50", "NA")
        default_phase_len = str(snakemake.params.default_phasing_min_length)
        default_core = str(snakemake.params.default_core_depth_threshold)

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
            "primary_mapped_reads": primary,
            "primary_mapped_percent": pct(primary, ns["reads"]),
            "mapq_ge_20_reads": q20,
            "mapq_ge_20_of_primary_percent": pct(q20, primary),
            "mapq_ge_30_reads": q30,
            "mapq_ge_30_of_primary_percent": pct(q30, primary),
            "mapq_ge_50_reads": q50,
            "mapq_ge_50_of_primary_percent": pct(q50, primary),
            "secondary_alignments": aqc.get("secondary", "NA"),
            "supplementary_alignments": aqc.get("supplementary", "NA"),
            "mean_length_mapq_ge_30": aqc.get("mean_length_mapq_ge_30", "NA"),
            "mean_length_mapq_lt_20": aqc.get("mean_length_mapq_lt_20", "NA"),
            "raw_alignment_mapped_percent": raw_mapped_pct,
            "raw_reference_coverage_percent": raw_cov_pct,
            "raw_mean_depth": raw_depth,
            "analysis_reference_coverage_percent": analysis_cov_pct,
            "analysis_mean_depth": analysis_depth,
            "observed_core_interval": core_region,
            "observed_core_length": core_length,
            "core_depth_threshold": meta.get("core_depth_threshold", "") or core_threshold or default_core,
            "target_region": meta.get("target_region", "") or "NONE",
            "phasing_min_length": meta.get("phasing_min_length", "") or default_phase_len,
        })
