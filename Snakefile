configfile: "config/config.yaml"

import csv
import re
from collections import defaultdict
from pathlib import Path

SAMPLES_TSV = config.get("samples", "config/samples.tsv")

SAFE_ID = re.compile(r"^[A-Za-z0-9._-]+$")
FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")


def clean_id(value, column):
    value = value.strip()
    if not value or not SAFE_ID.match(value):
        raise ValueError(
            f"Invalid {column}={value!r}. Use only letters, numbers, '.', '_' and '-'."
        )
    return value


rows = []
with open(SAMPLES_TSV, newline="") as handle:
    reader = csv.DictReader(handle, delimiter="\t")
    required = {"sample", "gene", "amplicon", "fastq_input", "reference"}
    missing = required - set(reader.fieldnames or [])
    if missing:
        raise ValueError(f"Missing columns in {SAMPLES_TSV}: {sorted(missing)}")

    for raw in reader:
        if not raw.get("sample", "").strip():
            continue
        row = {k: (v.strip() if isinstance(v, str) else v) for k, v in raw.items()}
        row["sample"] = clean_id(row["sample"], "sample")
        row["gene"] = clean_id(row["gene"], "gene")
        row["amplicon"] = clean_id(row["amplicon"], "amplicon")
        row["analysis"] = f'{row["sample"]}__{row["gene"]}'
        row["unit"] = f'{row["analysis"]}__{row["amplicon"]}'
        rows.append(row)

if not rows:
    raise ValueError(f"No amplicon rows found in {SAMPLES_TSV}")

UNITS = {}
ANALYSES = defaultdict(list)
ANALYSIS_META = {}

for row in rows:
    unit = row["unit"]
    analysis = row["analysis"]
    if unit in UNITS:
        raise ValueError(f"Duplicate sample/gene/amplicon combination: {unit}")
    UNITS[unit] = row
    ANALYSES[analysis].append(unit)

    previous = ANALYSIS_META.get(analysis)
    if previous and previous["reference"] != row["reference"]:
        raise ValueError(
            f"All amplicons in {analysis} must use the same reference. "
            f"Found {previous['reference']} and {row['reference']}"
        )
    ANALYSIS_META[analysis] = {
        "sample": row["sample"],
        "gene": row["gene"],
        "reference": row["reference"],
    }

UNIT_IDS = sorted(UNITS)
ANALYSIS_IDS = sorted(ANALYSES)
for analysis in ANALYSIS_IDS:
    ANALYSES[analysis] = sorted(ANALYSES[analysis])


def _is_fastq(path):
    s = str(path).lower()
    return any(s.endswith(ext) for ext in FASTQ_SUFFIXES)


def discover_fastqs(raw_input):
    p = Path(raw_input)
    if p.is_file():
        if not _is_fastq(p):
            raise ValueError(f"Input file is not FASTQ/FASTQ.GZ: {p}")
        return [str(p)]
    if p.is_dir():
        files = sorted(str(x) for x in p.rglob("*") if x.is_file() and _is_fastq(x))
        if not files:
            raise ValueError(f"No .fastq/.fastq.gz/.fq/.fq.gz files found under: {p}")
        return files
    raise FileNotFoundError(
        f"fastq_input does not exist: {p}. Check config/samples.tsv and your working directory."
    )


def unit_fastqs(wc):
    return discover_fastqs(UNITS[wc.unit]["fastq_input"])


def analysis_for_unit(unit):
    return UNITS[unit]["analysis"]


def analysis_reference(wc):
    return ANALYSIS_META[wc.analysis]["reference"]


def staged_reference_for_unit(wc):
    return f"results/reference/{analysis_for_unit(wc.unit)}/reference.fasta"


def staged_mmi_for_unit(wc):
    return f"results/reference/{analysis_for_unit(wc.unit)}/reference.mmi"


def reads_for_mapping(wc):
    if use_filtered:
        return f"results/filtered/{wc.unit}/{wc.unit}.filtered.fastq.gz"
    return f"results/input/{wc.unit}/{wc.unit}.combined.fastq.gz"


def amplicon_bams_for_analysis(wc):
    return [f"results/mapping/amplicons/{u}/{u}.sorted.bam" for u in ANALYSES[wc.analysis]]


use_filtered = bool(config.get("filtering", {}).get("enabled", False))
run_variants = bool(config.get("workflow", {}).get("run_variant_calling", False))
run_phasing = bool(config.get("workflow", {}).get("run_phasing", False)) and run_variants
run_consensus = bool(config.get("workflow", {}).get("run_consensus", False)) and run_phasing

include: "workflow/rules/input.smk"
include: "workflow/rules/qc.smk"
include: "workflow/rules/mapping.smk"
include: "workflow/rules/variants.smk"
include: "workflow/rules/phasing.smk"
include: "workflow/rules/consensus.smk"
include: "workflow/rules/report.smk"

final_targets = ["results/summary/input_manifest.tsv"]
final_targets += expand("results/input/{unit}/{unit}.combined.fastq.gz", unit=UNIT_IDS)
final_targets += expand("results/qc/raw/{unit}/NanoStats.txt", unit=UNIT_IDS)
final_targets += expand("results/mapping/amplicons/{unit}/{unit}.sorted.bam.bai", unit=UNIT_IDS)
final_targets += expand("results/mapping/amplicons/{unit}/{unit}.coverage.txt", unit=UNIT_IDS)
final_targets += expand("results/mapping/genes/{analysis}/{analysis}.merged.bam.bai", analysis=ANALYSIS_IDS)
final_targets += expand("results/mapping/genes/{analysis}/{analysis}.coverage.txt", analysis=ANALYSIS_IDS)
final_targets += [
    "results/summary/amplicon_summary.tsv",
    "results/summary/gene_summary.tsv",
]

if run_variants:
    final_targets += expand("results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi", analysis=ANALYSIS_IDS)
if run_phasing:
    final_targets += expand("results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi", analysis=ANALYSIS_IDS)
if run_consensus:
    final_targets += expand("results/consensus/{analysis}/{analysis}.haplotype1.fasta", analysis=ANALYSIS_IDS)
    final_targets += expand("results/consensus/{analysis}/{analysis}.haplotype2.fasta", analysis=ANALYSIS_IDS)


rule all:
    input:
        final_targets
