"""Manifest identity and transparent technical-QC policy for Excel reporting.

This module never changes workflow analysis keys or production thresholds.
Reporting groups are (run ID, workflow sample ID, display sample ID, barcode); the existing analysis
and unit keys remain ``sample__gene`` and ``sample__gene__amplicon``.  If those
keys pool distinct reporting identities, callers must withhold pooled metrics.

Policy schema (no cutoffs are installed by this module)::

    sample_types:
      sample:
        criteria:
          mean_target_depth:
            min: 100  # example only; not a validated default
            required: true
            severity: REVIEW
      negative_control:
        criteria: {}  # must be configured separately
    genes:
      ABO:
        sample_types:
          sample:
            criteria:
              phased_fraction:
                min: 1
                not_applicable_if: {eligible_heterozygous_records: 0}

Gene-specific criteria override same-named sample-type criteria.  A null
override removes that criterion. Controls never inherit ordinary sample rules.
Only min/max/equals/one_of comparisons are supported; missing thresholds are
NOT ASSESSED, and an actual numerical zero remains a measured value.
"""

from __future__ import annotations

import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path


STATUSES = ("PASS", "REVIEW", "NOT ASSESSED", "FAIL")
PRECEDENCE = {status: index for index, status in enumerate(STATUSES)}
SAMPLE_TYPES = {"sample", "positive_control", "negative_control"}
OPERATORS = {"min", "max", "equals", "one_of"}


def load_manifest(path, run_id=""):
    """Return every manifest-defined reporting group with original text IDs.

    Optional columns are sample_id (display label), barcode, run_id and
    sample_type.  The legacy sample column remains the analysis identifier.
    Missing barcode/run IDs are empty strings with explicit metadata notes;
    neither is inferred from paths or FASTQ headers. Incompatible metadata and
    exact duplicate rows are errors. Reused workflow keys across identities
    remain separate reporting groups and are explicitly flagged as ambiguous.
    """
    groups = {}
    analysis_groups = defaultdict(set)
    unit_identities = {}
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "gene", "amplicon", "fastq_input", "reference"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing manifest columns: {', '.join(sorted(missing))}")
        if len(reader.fieldnames) != len(set(reader.fieldnames)):
            raise ValueError("Duplicate column names in sample manifest")
        for line, raw in enumerate(reader, start=2):
            if None in raw:
                raise ValueError(f"Manifest line {line} has extra tab-separated fields")
            if not (raw.get("sample") or "").strip():
                continue  # Match the production Snakefile's blank-sample behavior.
            row = {key: (value or "").strip() for key, value in raw.items()}
            for column in required:
                if not row[column]:
                    raise ValueError(f"Manifest line {line}: empty {column}")
            for column in ("sample", "gene", "amplicon"):
                if not re.fullmatch(r"[A-Za-z0-9._-]+", row[column]):
                    raise ValueError(f"Manifest line {line}: invalid workflow {column}; use sample_id for display labels")
            sample_id = row.get("sample_id") or row["sample"]
            barcode = row.get("barcode", "")
            row_run_id = row.get("run_id") or str(run_id or "")
            sample_type = row.get("sample_type") or "sample"
            if sample_type not in SAMPLE_TYPES:
                raise ValueError(
                    f"Manifest line {line}: unknown sample_type {sample_type!r}; "
                    f"expected one of {', '.join(sorted(SAMPLE_TYPES))}"
                )
            row.update(
                sample_id=sample_id, barcode=barcode, run_id=row_run_id,
                sample_type=sample_type, manifest_line=line,
                analysis=f"{row['sample']}__{row['gene']}",
                unit=f"{row['sample']}__{row['gene']}__{row['amplicon']}",
            )
            # Display labels never group distinct workflow samples, including
            # repeated Positive/Negative labels when no barcode is supplied.
            key = (row_run_id, row["sample"], sample_id, barcode)
            if (row["unit"], key) in unit_identities:
                raise ValueError(
                    f"Duplicate workflow unit {row['unit']!r} at manifest lines "
                    f"{unit_identities[row['unit'], key]} and {line} for one reporting identity"
                )
            unit_identities[row["unit"], key] = line
            if key not in groups:
                groups[key] = {
                    "sample_key": json.dumps(key, ensure_ascii=False, separators=(",", ":")),
                    "run_id": row_run_id, "sample_id": sample_id, "barcode": barcode,
                    "workflow_sample": row["sample"],
                    "sample_type": sample_type, "genes": [], "rows": [], "issues": [],
                    "identity_notes": [], "ambiguous_analyses": [],
                }
                if not barcode:
                    groups[key]["identity_notes"].append("Barcode not provided in manifest")
                if not row_run_id:
                    groups[key]["identity_notes"].append("Run ID not provided")
            group = groups[key]
            if group["sample_type"] != sample_type:
                raise ValueError(
                    f"Conflicting sample_type values for reporting identity {key!r}"
                )
            if row["gene"] not in group["genes"]:
                group["genes"].append(row["gene"])
            group["rows"].append(row)
            analysis_groups[row["analysis"]].add(key)
    if not groups:
        raise ValueError(f"No sample rows found in {path}")
    for analysis, keys in analysis_groups.items():
        if len(keys) > 1:
            reason = (
                f"Analysis {analysis} pools {len(keys)} distinct manifest sample/barcode/run "
                "identities; gene-level outputs cannot be attributed to one reporting row"
            )
            for key in keys:
                groups[key]["issues"].append(reason)
                groups[key]["ambiguous_analyses"].append(analysis)
    for group in groups.values():
        by_gene = defaultdict(set)
        for row in group["rows"]:
            by_gene[row["gene"]].add(row["analysis"])
        for gene, analyses in by_gene.items():
            if len(analyses) > 1:
                group["issues"].append(
                    f"Reporting identity has multiple analyses for {gene}: "
                    + ", ".join(sorted(analyses))
                    + "; results cannot be silently pooled"
                )
                group["ambiguous_analyses"].extend(sorted(analyses))
        group["ambiguous_analyses"] = sorted(set(group["ambiguous_analyses"]))
        for row in group["rows"]:
            row["analysis_identity_ambiguous"] = row["analysis"] in group["ambiguous_analyses"]
    return list(groups.values())


def resolve_criteria(policy, sample_type="sample", gene=""):
    """Resolve criteria using only the specified sample/control type."""
    if sample_type not in SAMPLE_TYPES:
        raise ValueError(f"Unknown sample_type {sample_type!r}")
    policy = policy or {}
    if not isinstance(policy, dict):
        raise ValueError("QC policy must be a mapping")
    criteria = {}
    locations = [policy.get("sample_types", {}).get(sample_type, {})]
    locations.append(policy.get("genes", {}).get(gene, {}).get("sample_types", {}).get(sample_type, {}))
    for location in locations:
        if not isinstance(location, dict):
            raise ValueError(f"Policy for {sample_type}/{gene} must be a mapping")
        configured = location.get("criteria", {})
        if not isinstance(configured, dict):
            raise ValueError(f"criteria for {sample_type}/{gene} must be a mapping")
        for metric, rule in configured.items():
            if rule is None:
                criteria.pop(metric, None)
            elif not isinstance(rule, dict):
                raise ValueError(f"Criterion {metric} must be a mapping or null")
            else:
                criteria[metric] = dict(rule)
    return criteria


def _is_missing(value):
    return value is None or (isinstance(value, float) and not math.isfinite(value))


def _condition_matches(metrics, condition):
    if not isinstance(condition, dict) or not condition:
        raise ValueError("not_applicable_if must be a nonempty metric:value mapping")
    return all(key in metrics and not _is_missing(metrics[key]) and metrics[key] == value
               for key, value in condition.items())


def evaluate_criterion(metric, value, rule, metrics=None):
    """Evaluate one declared metric and keep its evidence and threshold visible.

    required controls whether absent evidence/criteria prevent PASS. A failed
    optional criterion is REVIEW even if its configured severity is FAIL.
    """
    allowed = OPERATORS | {"required", "severity", "not_applicable_if", "description", "units"}
    unknown = set(rule) - allowed
    if unknown:
        raise ValueError(f"Unknown criterion options for {metric}: {sorted(unknown)}")
    required = rule.get("required", True)
    if not isinstance(required, bool):
        raise ValueError(f"Criterion {metric}: required must be a boolean")
    severity = rule.get("severity", "FAIL")
    if severity not in {"FAIL", "REVIEW"}:
        raise ValueError(f"Criterion {metric}: severity must be FAIL or REVIEW")
    for operator in ("min", "max"):
        if rule.get(operator) is not None and (isinstance(rule[operator], bool)
                                 or not isinstance(rule[operator], (int, float))
                                 or not math.isfinite(rule[operator])):
            raise ValueError(f"Criterion {metric}: {operator} must be a finite number")
    if rule.get("min") is not None and rule.get("max") is not None and rule["min"] > rule["max"]:
        raise ValueError(f"Criterion {metric}: min is greater than max")
    if rule.get("one_of") is not None and not isinstance(rule["one_of"], (list, tuple)):
        raise ValueError(f"Criterion {metric}: one_of must be a list")
    result = {"metric": metric, "value": None if _is_missing(value) else value,
              "required": required, "rule": dict(rule)}
    if metric in (metrics or {}).get("_not_applicable", []):
        result.update(status="NOT APPLICABLE", reason=f"{metric}: not applicable to the observed data")
        return result
    if "not_applicable_if" in rule and _condition_matches(metrics or {}, rule["not_applicable_if"]):
        result.update(status="NOT APPLICABLE", reason=f"{metric}: not applicable ({rule['not_applicable_if']})")
        return result
    comparisons = {operator: rule[operator] for operator in sorted(OPERATORS) if rule.get(operator) is not None}
    if not comparisons:
        result.update(status="NOT ASSESSED", reason=f"{metric}: no evaluation threshold or expected value configured")
        return result
    if _is_missing(value):
        result.update(status="NOT ASSESSED", reason=f"{metric}: required evidence is unavailable" if required
                      else f"{metric}: optional evidence is unavailable")
        return result
    failures = []
    if "min" in comparisons or "max" in comparisons:
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"Criterion {metric}: numeric comparison requires a number, got {value!r}")
    if "min" in comparisons and value < rule["min"]:
        failures.append(f"{value} < minimum {rule['min']}")
    if "max" in comparisons and value > rule["max"]:
        failures.append(f"{value} > maximum {rule['max']}")
    if "equals" in comparisons and value != rule["equals"]:
        failures.append(f"{value!r} != expected {rule['equals']!r}")
    if "one_of" in comparisons and value not in rule["one_of"]:
        failures.append(f"{value!r} is outside {rule['one_of']!r}")
    if failures:
        result.update(status=severity if required else "REVIEW",
                      reason=f"{metric}: " + "; ".join(failures))
    else:
        result.update(status="PASS", reason=f"{metric}: {value!r} satisfies {comparisons}")
    return result


def aggregate_assessments(assessments):
    """Aggregate gene/criterion assessments without concealing any explanation.

    FAIL > NOT ASSESSED > REVIEW > PASS. NOT APPLICABLE is excluded. Optional
    missing criteria do not prevent PASS but their explanations remain visible.
    A wholly inapplicable/unconfigured assessment is never promoted to PASS.
    """
    reasons, statuses = [], []
    for item in assessments:
        status = item["status"]
        if status not in PRECEDENCE and status != "NOT APPLICABLE":
            raise ValueError(f"Unknown assessment status {status!r}")
        for reason in item.get("reasons", [item.get("reason", "")]):
            if reason and reason not in reasons:
                reasons.append(reason)
        if status == "NOT APPLICABLE":
            continue
        if status == "NOT ASSESSED" and not item.get("required", True):
            continue
        statuses.append(status)
    if statuses:
        status = max(statuses, key=PRECEDENCE.__getitem__)
    else:
        status = "NOT ASSESSED"
        reasons.append("No applicable criteria were assessed")
    return {"status": status, "reasons": reasons}


def evaluate_gene(metrics, policy, sample_type="sample", gene="", evidence_issues=None):
    """Return status, all reasons and every criterion for one expected gene.

    evidence_issues are independent required-evidence failures from collection
    (missing files, parse errors, ambiguous identity); each prevents PASS without
    converting missing measurements to zero. They do not imply biological FAIL.
    """
    criteria = resolve_criteria(policy, sample_type, gene)
    assessed = [evaluate_criterion(metric, metrics.get(metric), rule, metrics)
                for metric, rule in criteria.items()]
    if not criteria:
        assessed.append({"metric": "evaluation_policy", "value": None, "required": True,
                         "status": "NOT ASSESSED", "rule": {},
                         "reason": f"No evaluation criteria configured for {sample_type}/{gene}"})
    for issue in evidence_issues or []:
        assessed.append({"metric": "required_evidence", "value": None, "required": True,
                         "status": "NOT ASSESSED", "rule": {}, "reason": str(issue)})
    result = aggregate_assessments(assessed)
    result["criteria"] = assessed
    return result
