#!/usr/bin/env python3
"""Final Snakemake evaluation workbook, or explicit standalone partial report."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys

if "snakemake" in globals():
    sys.path.insert(0, str(snakemake.params.script_dir))
else:
    sys.path.insert(0, str(Path(__file__).resolve().parent))

import yaml
from evaluation_metrics import Evidence, collect_gene, table, workflow_stages
from evaluation_policy import load_manifest, evaluate_gene, aggregate_assessments
from evaluation_workbook import write_workbook

IDENTITY = ["run_id", "sample_id", "barcode", "sample_type", "gene", "analysis", "workflow_sample"]
HEADERS = {
    "Sample_Summary": ["run_id", "sample_id", "barcode", "sample_type", "expected_genes", "execution_status", "overall_qc", "gene_results", "reasons"],
    "Gene_QC": IDENTITY + ["execution_status", "gene_qc", "reasons", "input_reads", "mapper_input_reads", "primary_mapped_reads", "mapping_percent", "retained_variant_reads", "phasing_reads", "target_bases", "mean_target_depth", "median_target_depth", "min_target_depth", "assessment_depth", "target_breadth_fraction", "zero_coverage_bases", "zero_coverage_regions", "low_coverage_bases", "low_coverage_regions", "expected_amplicons", "observed_amplicons", "variant_records", "multiallelic_records", "eligible_heterozygous_records", "phased_heterozygous_records", "unphased_heterozygous_records", "phased_fraction", "phase_block_count", "assigned_read_fraction", "accepted_variants", "unresolved_variants", "hp1_target_unmasked_fraction", "hp2_target_unmasked_fraction"],
    "Variants": IDENTITY + ["contig", "position", "REF", "ALT", "genotype", "FILTER", "quality", "depth", "allele_depths", "allele_fractions", "phase_set"],
    "Haplotypes": IDENTITY + ["row_type", "haplotype", "contig", "phase_set", "consensus_length", "target_masked_reference_bases", "target_unmasked_fraction"],
    "Run_Info": ["category", "name", "value", "definition", "source_file"],
}


def merge_config(target, overlay):
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(target.get(key), dict):
            merge_config(target[key], value)
        else:
            target[key] = value
    return target


def identity(group):
    return {k: group[k] or None for k in ("run_id", "sample_id", "barcode", "sample_type", "workflow_sample")}


def execution_status(sources, issues, strict):
    if any(s["state"] == "PARSE ERROR" for s in sources) or any(i["status"] == "FAIL" for i in issues):
        return "EVIDENCE ERROR (analysis execution unverified)"
    if any(s["state"] == "MISSING" and s["required"] for s in sources):
        return "INCOMPLETE: required outputs missing"
    return "DEPENDENCIES SATISFIED (cached outputs may be reused)" if strict else "OUTPUTS AVAILABLE; EXECUTION NOT VERIFIED"


def combine_execution(values, strict):
    if any("ERROR" in v for v in values):
        return "EVIDENCE ERROR (see gene details)"
    if any("INCOMPLETE" in v for v in values):
        return "INCOMPLETE"
    if any("UNVERIFIED" in v or "NOT VERIFIED" in v for v in values):
        return "PARTIAL REPORT: EXECUTION NOT VERIFIED"
    return "DEPENDENCIES SATISFIED (cached outputs may be reused)" if strict else "PARTIAL REPORT: EXECUTION NOT VERIFIED"


def build_report(config, workdir, output, repository, *, strict=False, config_files=(), dependencies=(), synthetic=False):
    """Generate a workbook/JSON; strict=True is reserved for Snakemake context.

    Standalone presence of old files never establishes successful execution.
    A partial report can assess readable historical metrics, while its overall
    QC remains NOT ASSESSED unless a configured FAIL already takes precedence.
    """
    workdir, repository, output = Path(workdir).resolve(), Path(repository).resolve(), Path(output).resolve()
    manifest = Path(config.get("samples", "config/samples.tsv"))
    manifest = manifest if manifest.is_absolute() else workdir / manifest
    groups = load_manifest(manifest, config.get("reporting", {}).get("run_id") or "")
    results = workdir / "results"
    shared = Evidence()
    input_manifest_path = results / "summary/input_manifest.tsv"
    input_manifest = shared.read(input_manifest_path, lambda p: table(p, ["sample", "gene", "amplicon", "unit", "source"]), stage="input provenance")
    # Every strict dependency was selected by the current Snakemake DAG. This
    # confirms completion of prerequisites, not that cached outputs ran anew.
    missing_deps = [str(p) for p in dependencies if not Path(p).exists()]
    if strict and missing_deps:
        raise ValueError("Missing strict report dependencies: " + "; ".join(missing_deps))
    sheets = {name: [] for name in HEADERS}
    details, sources, strict_errors = [], list(shared.sources), []
    policy = config.get("reporting", {}).get("qc", {})
    for group in groups:
        gene_assessments, gene_executions, results_text = [], [], []
        for gene in group["genes"]:
            rows = [r for r in group["rows"] if r["gene"] == gene]
            metrics, variants, haplotypes, evidence = collect_gene(rows, results, config, input_manifest)
            issues = shared.issues + evidence.issues
            issue_reasons = [f"{i['status']}: {i['reason']}" for i in issues] + group["issues"]
            unrequired_control_evidence = []
            if group["sample_type"] in {"positive_control", "negative_control"}:
                # Control requirements come from their own configured criteria.
                # Missing ordinary-sample outputs remain execution gaps, without
                # silently becoming coverage/phasing requirements for controls.
                issue_reasons = [f"{i['status']}: {i['reason']}" for i in issues if i["status"] == "FAIL"] + group["issues"]
                unrequired_control_evidence = [i["reason"] + " (output unavailable; control QC uses only its configured criteria)" for i in issues if i["status"] != "FAIL"]
            state = execution_status(shared.sources + evidence.sources, issues, strict)
            if any(r.get("analysis_identity_ambiguous") for r in rows):
                state = "AMBIGUOUS IDENTITY: EXECUTION NOT VERIFIED"
            assessment = evaluate_gene(metrics, policy, sample_type=group["sample_type"], gene=gene, evidence_issues=issue_reasons)
            assessment["reasons"].extend(unrequired_control_evidence)
            gene_assessments.append({**assessment, "reasons": [f"{gene}: {r}" for r in assessment["reasons"]]})
            gene_executions.append(state)
            values = {k: metrics.get(k) if metrics.get(k) is not None else 'NOT AVAILABLE' for k in ('variant_records', 'accepted_variants', 'unresolved_variants')}
            results_text.append(f"{gene}: {assessment['status']}; calls={values['variant_records']}; accepted={values['accepted_variants']}; unresolved={values['unresolved_variants']}")
            row = {**identity(group), "gene": gene, "analysis": metrics["analysis"], "execution_status": state,
                   "gene_qc": assessment["status"], "reasons": " | ".join(assessment["reasons"])}
            row.update({k: v for k, v in metrics.items() if not k.startswith("_")})
            for key in metrics.get("_not_applicable", []):
                if row.get(key) is None:
                    row[key] = "NOT APPLICABLE"
            row["metrics_file"] = str(output.with_suffix(".json"))
            sheets["Gene_QC"].append(row)
            for record in variants:
                sheets["Variants"].append({**identity(group), "gene": gene, "analysis": metrics["analysis"], **record})
            for record in haplotypes:
                sheets["Haplotypes"].append({**identity(group), "gene": gene, "analysis": metrics["analysis"], **record})
            if not haplotypes:
                sheets["Haplotypes"].append({**identity(group), "gene": gene, "analysis": metrics["analysis"],
                                            "row_type": "AVAILABILITY", "applicability": metrics.get("phasing_applicability", "NOT ASSESSED"),
                                            "reasons": "No haplotype/block evidence available; FASTA existence alone does not establish phase"})
            sources.extend(evidence.sources)
            strict_errors.extend(s for s in evidence.sources + shared.sources if s["required"] and s["state"] != "READ")
            strict_errors.extend(i for i in issues if i["status"] == "FAIL")
            details.append({"identity": identity(group), "gene": gene, "metrics": metrics,
                            "assessment": assessment, "execution_status": state, "issues": issues, "sources": evidence.sources})
        aggregation_inputs = list(gene_assessments)
        if not strict:
            aggregation_inputs.append({"status": "NOT ASSESSED", "reasons": ["Partial report: current execution completion is not verified; available files may be stale"], "criteria": []})
        combined = aggregate_assessments(aggregation_inputs)
        summary = {**identity(group), "expected_genes": ", ".join(group["genes"]),
                   "execution_status": combine_execution(gene_executions, strict), "overall_qc": combined["status"],
                   "gene_results": " | ".join(results_text), "reasons": " | ".join(combined["reasons"] + group["identity_notes"]),
                   "report_mode": "SYNTHETIC DEMONSTRATION" if synthetic else "STRICT DEPENDENCY REPORT" if strict else "PARTIAL / INCOMPLETE",
                   "manifest_file": str(manifest), "evidence_file": str(output.with_suffix(".json"))}
        for gene, assessment in zip(group["genes"], gene_assessments):
            summary[f"{gene}_QC"] = assessment["status"]
        sheets["Sample_Summary"].append(summary)
    if strict_errors and strict:
        # Do not leave a normal-looking workbook after dependency or parse errors.
        raise ValueError("Strict report rejected unreadable/inconsistent evidence: " + json.dumps(strict_errors))
    info = sheets["Run_Info"]
    def note(category, name, value, definition="", source_file=None):
        if isinstance(value, (dict, list)):
            value = json.dumps(value, ensure_ascii=False, sort_keys=True)
        info.append({"category": category, "name": name, "value": value, "definition": definition, "source_file": source_file})
    created = datetime.now(timezone.utc).isoformat()
    try:
        commit = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=repository, text=True, stderr=subprocess.DEVNULL).strip()
        dirty = bool(subprocess.check_output(["git", "status", "--porcelain"], cwd=repository, text=True))
    except (OSError, subprocess.CalledProcessError):
        commit, dirty = None, None
    note("Run", "report_mode", "SYNTHETIC DEMONSTRATION" if synthetic else "STRICT DEPENDENCY REPORT" if strict else "PARTIAL / INCOMPLETE")
    note("Run", "generated_utc", created)
    note("Run", "work_directory", str(workdir), "Analysis output namespace; not inferred sequencing-run identity")
    note("Run", "configured_run_id", config.get("reporting", {}).get("run_id"), "Manifest run_id overrides this reporting fallback")
    note("Provenance", "report_code_git_commit", commit, "Report-generation checkout; not proof historical outputs used this commit")
    note("Provenance", "report_checkout_has_uncommitted_changes", dirty)
    note("Provenance", "original_analysis_git_commit", None, "Not recorded in the standard analysis artifacts; do not substitute report-generation commit")
    note("Provenance", "manifest_sha256", hashlib.sha256(manifest.read_bytes()).hexdigest(), source_file=str(manifest))
    note("Tools", "python", platform.python_version(), "Report-generation interpreter only")
    note("Tools", "architecture", platform.machine(), "Report-generation host architecture")
    for name in ("evaluation_report.py", "evaluation_metrics.py", "evaluation_policy.py", "evaluation_workbook.py"):
        path = repository / "workflow/scripts" / name
        note("Provenance", "report_source_sha256", hashlib.sha256(path.read_bytes()).hexdigest(), source_file=str(path))
    for p in config_files:
        path = Path(p).resolve()
        note("Provenance", "config_sha256", hashlib.sha256(path.read_bytes()).hexdigest(), source_file=str(path))
    for name in ("openpyxl", "PyYAML", "snakemake", "pysam", "whatshap"):
        try:
            version = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            version = None
        note("Tools", name, version, "Report-generation environment only; original analysis tool version not inferred")
    for key in ("workflow", "filtering", "mapping", "clair3", "phasing", "haplotypes", "reporting"):
        note("Configuration", key, config.get(key), "Current effective requested configuration; partial reports cannot attest historical parameters")
    definitions = {
        "PASS": "All required configured gene criteria assessed and passed; technical QC only, never biological validation.",
        "REVIEW": "Configured warning condition needs review, with no higher-precedence condition.",
        "FAIL": "A required configured criterion failed; zero variants alone never fails.",
        "NOT ASSESSED": "Required evidence/evaluation criteria unavailable; missing required evidence prevents PASS.",
        "aggregation": "FAIL > NOT ASSESSED > REVIEW > PASS. Every expected gene participates; all reasons retained. N/A criteria excluded. Partial overall status cannot PASS.",
        "NOT APPLICABLE": "Criterion does not apply, e.g. phase fraction without eligible heterozygotes; distinct from unavailable evidence.",
        "controls": "sample, positive_control and negative_control policies are separate; controls never inherit ordinary-sample coverage requirements. Absent control policy is NOT ASSESSED.",
        "identifiers": "Manifest text preserved. Optional sample_id is a display ID; sample remains the unique workflow ID. No barcodes inferred. Ambiguous pooled analysis keys withheld.",
        "missing_values": "NOT AVAILABLE cells are unavailable annotations/metrics; actual zero remains numeric zero. Parse errors are explicit evidence errors, not ordinary missing values.",
        "input_reads": "FASTQ observations from raw NanoStats, summed once across distinct source files in the gene; no totals pooled across genes or units with repeated source paths.",
        "mapping_percent": "100 * primary_mapped / primary_records from alignment_qc; primary_records includes unmapped mapper inputs and excludes secondary/supplementary records. May differ from original FASTQ count when filtering enabled.",
        "retained_reads": "Gene variant/phasing flagstat primary mapped observations. BAMs retain alignment outside target; primary/MAPQ filters apply, phasing additionally uses configured per-amplicon minimum length. Not independent molecules.",
        "observed_amplicons": "Number of expected manifest amplicons with readable alignment QC, not proof of amplicon identity or a minimum-yield pass. Per-amplicon input/mapped observations are preserved in the JSON evidence sidecar.",
        "target_depth": "Per-amplicon samtools depth -aa arrays summed over merged manifest target union; all positions including zeros required. Aligned bases only, no deletions/skips; default depth excludes unmapped/secondary/DUP/QCFAIL; filtered inputs also exclude supplementary. BQ minimum 0; no pileup cap.",
        "target_breadth_fraction": "Fraction of union reference positions >= assessment_depth (default existing callable_min_depth). Measurement threshold is not itself a validated gene acceptance fraction.",
        "coordinates": "VCF positions and manifest targets 1-based; manifest ends inclusive. BED regions in JSON are 0-based half-open. Coordinates refer to staged reference identifiers, not assumed genome assembly.",
        "multiallelic": "One Variants row per original normalized record. ALT order retained; AD is REF then ALT vector; AF follows ALT order. Missing vector elements remain '.'; DP need not equal sum(AD). No implicit splitting or AF imputation.",
        "phasing": "Eligible heterozygotes are biallelic phased-VCF records. PS plus phased GT required. Phase blocks local to contig/gene; HP labels do not establish linkage between disconnected blocks or genes. No heterozygotes is N/A, not failure.",
        "consensus_completeness": "Unmasked intended reference positions / intended reference positions, from BED masks. FASTA length/N counts are separate; indels shift sequence coordinates. Supported deletions can be resolved without aligned-base depth. This metric is not phase or allele accuracy.",
        "support_decisions": "ACCEPT/UNRESOLVED are existing workflow decisions, not automatically gene QC or proof of variant truth. Accepted alleles can still be masked for low depth.",
        "limitations": "No phenotype or named-allele interpretation implemented. Haplotype imbalance alone is not allele dropout. No truth-set or assay validation inferred. Historical outputs/logs do not prove successful current execution.",
        "unavailable_metrics": "No unique-molecule counts, biological genotype truth, validated gene-wide thresholds, clinical interpretation, or original tool versions unless explicit provenance is supplied. HP median depth is not present in existing structured HP QC and is not invented.",
    }
    for name, description in definitions.items():
        note("Definitions", name, description)
    for item in details:
        note("Reference", item["metrics"]["analysis"], item["metrics"].get("reference_identifiers"),
             "Identifiers in the staged FASTA; coordinate system for this gene", item["metrics"].get("reference_file"))
        note("Targets", item["metrics"]["analysis"], item["metrics"].get("target_intervals"),
             "Merged intended manifest union; 1-based inclusive coordinates")
        for artifact in item["metrics"].get("_vcf_provenance", []):
            for header in artifact["headers"] or []:
                note("Recorded analysis metadata", item["metrics"]["analysis"], header, "Literal VCF header provenance; distinct from report-environment versions", artifact["file"])
        for criterion in item["assessment"]["criteria"]:
            label = (f"{item['metrics']['analysis']} / {item['identity']['sample_id']} / "
                     f"run={item['identity']['run_id'] or 'NOT AVAILABLE'} / "
                     f"barcode={item['identity']['barcode'] or 'NOT AVAILABLE'}")
            note("Assessment", label, criterion)
    seen = set()
    for source in sources:
        key = source["file"], source["stage"]
        if key not in seen:
            note("Metric sources", source["stage"], {k:v for k,v in source.items() if k != "file"}, source_file=source["file"])
            seen.add(key)
    payload = {"schema_version": 1, "generated_utc": created, "mode": "strict" if strict else "partial", "synthetic": synthetic,
               "manifest": str(manifest), "workdir": str(workdir), "config": config, "samples": sheets["Sample_Summary"],
               "genes": details, "sources": sources, "definitions": definitions}
    output.parent.mkdir(parents=True, exist_ok=True)
    sidecar = output.with_suffix(".json")
    temporary = sidecar.with_suffix(".json.tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    temporary.replace(sidecar)
    write_workbook(output, sheets, metadata={"headers": HEADERS, "link_base_dir": str(workdir),
                                            "title": "SYNTHETIC DEMONSTRATION" if synthetic else "Pipeline evaluation"})
    return payload


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--partial", action="store_true", required=True, help="Read available outputs; never assert current execution success")
    p.add_argument("--config", nargs="+", type=Path, required=True, help="YAML configuration overlays in application order")
    p.add_argument("--workdir", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path)
    p.add_argument("--synthetic", action="store_true", help="Prominently label synthetic demonstrations")
    args = p.parse_args()
    repo = Path(__file__).resolve().parents[2]
    cfg = {}
    for path in args.config:
        document = yaml.safe_load(path.read_text())
        if not isinstance(document, dict):
            raise ValueError(f"Configuration must be a mapping: {path}")
        merge_config(cfg, document)
    output = args.output or args.workdir / "results/reports/pipeline_evaluation.partial.xlsx"
    report = build_report(cfg, args.workdir, output, repo, config_files=args.config, synthetic=args.synthetic)
    print(f"Partial workbook: {output}; expected sample/barcodes={len(report['samples'])}, expected sample genes={len(report['genes'])}")


if "snakemake" in globals():
    build_report(dict(snakemake.config), Path.cwd(), snakemake.output.xlsx, snakemake.params.repository,
                 strict=True, dependencies=list(snakemake.input))
elif __name__ == "__main__":
    main()
