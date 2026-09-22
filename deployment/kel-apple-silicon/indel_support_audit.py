#!/usr/bin/env python3
"""Historical read-only audit of the archived defective exact-indel counter.

Only the observations/count_support function definitions are executed. The
Snakemake script's top-level workflow and subprocess calls never execute.
JSON goes to stdout; callers may redirect it to a log. No repository file changes.
"""

import argparse
import ast
import hashlib
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--source", type=Path,
        default=Path(__file__).resolve().parent / "indel_fix/before/haplotype_variant_support.py",
    )
    args = parser.parse_args()
    source = args.source.read_text()
    tree = ast.parse(source, filename=str(args.source))
    names = {"observations", "count_support"}
    functions = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in names]
    if {node.name for node in functions} != names:
        raise RuntimeError("Required audit functions are missing from the source")
    namespace = {}
    exec(compile(ast.Module(body=functions, type_ignores=[]), str(args.source), "exec"), namespace)

    # Non-observations may alternatively be excluded from denominators, but they
    # must never be counted as exact REF/ALT. Here they are reported as OTHER.
    cases = [
        {
            "name": "deletion_placeholders_and_reference_skips_are_not_exact_ref",
            "ref": "AC", "alt": "A", "pileup_bases": "*#<>",
            "expected": {"ref": 0, "alt": 0, "other": 4, "supported": True},
        },
        {
            "name": "ambiguous_or_mismatching_indel_anchors_are_not_exact_ref",
            "ref": "AC", "alt": "A", "pileup_bases": "NnGg",
            "expected": {"ref": 0, "alt": 0, "other": 4, "supported": True},
        },
        {
            "name": "matching_insertion_with_wrong_anchor_is_not_exact_alt",
            "ref": "A", "alt": "AT", "pileup_bases": "G+1t",
            "expected": {"ref": 0, "alt": 0, "other": 1, "supported": True},
        },
        {
            "name": "matching_deletion_with_wrong_anchor_is_not_exact_alt",
            "ref": "AC", "alt": "A", "pileup_bases": "G-1c",
            "expected": {"ref": 0, "alt": 0, "other": 1, "supported": True},
        },
        {
            "name": "control_matching_anchor_and_insertion",
            "ref": "A", "alt": "AT", "pileup_bases": ".+1t",
            "expected": {"ref": 0, "alt": 1, "other": 0, "supported": True},
        },
        {
            "name": "control_matching_anchor_without_indel",
            "ref": "AC", "alt": "A", "pileup_bases": ".,",
            "expected": {"ref": 2, "alt": 0, "other": 0, "supported": True},
        },
    ]
    for case in cases:
        observed = namespace["count_support"](case["ref"], case["alt"], case["pileup_bases"])
        case["observed"] = dict(zip(["ref", "alt", "other", "supported"], observed))
        case["matches_expected"] = case["observed"] == case["expected"]
    print(json.dumps({
        "source": str(args.source.resolve()),
        "source_sha256": hashlib.sha256(args.source.read_bytes()).hexdigest(),
        "scope": "Function-level read-only audit; no workflow or biological validation performed",
        "interpretation": "Unexpected counts reproduce an exact-indel support classification defect",
        "cases": cases,
        "matching_cases": sum(case["matches_expected"] for case in cases),
        "mismatching_cases": sum(not case["matches_expected"] for case in cases),
    }, indent=2))


if __name__ == "__main__":
    main()
