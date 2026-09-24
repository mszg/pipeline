#!/usr/bin/env python3
"""Verify local KEL core outputs and source preservation, without changing reads."""
import csv
import hashlib
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
SETUP = ROOT / "setup"
RESULTS = ROOT / "workflow/results"
ANALYSIS = "KEL__KEL"


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    audit = json.loads((SETUP / "input_audit/kel_input_audit.json").read_text())
    report = {"scope": "KEL core outputs and original-file preservation", "sources": [], "amplicons": []}
    for entry in audit["files"]:
        path = Path(audit["source_directory"]) / entry["relative_path"]
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        stat = path.stat()
        assert digest == entry["sha256"], path
        assert stat.st_size == entry["compressed_bytes"], path
        assert stat.st_mtime_ns == entry["mtime_ns"], path
        report["sources"].append({"file": entry["relative_path"], "sha256_unchanged": True})

    manifest = table(RESULTS / "summary/input_manifest.tsv")
    assert len(manifest) == len(audit["files"])
    expected_sources = {
        (str((Path(audit["source_directory"]) / entry["relative_path"]).resolve()),
         "fragment" + Path(entry["relative_path"]).parts[0].rsplit("_", 1)[-1])
        for entry in audit["files"]
    }
    observed_sources = {
        (str((ROOT / "workflow" / row["source"]).resolve()), row["amplicon"])
        for row in manifest
    }
    assert observed_sources == expected_sources
    assert all(row["sample"] == row["gene"] == "KEL" and row["unit"] == f'{ANALYSIS}__{row["amplicon"]}' for row in manifest)
    target = RESULTS / f"targets/{ANALYSIS}/{ANALYSIS}.bed"
    assert target.read_text().strip() == "NG_007492.3\t410\t28023"
    report["target_union_bases"] = 27613
    report["overlap_bases"] = 1224
    summaries = table(RESULTS / "summary/amplicon_summary.tsv")
    aggregate_depth = [0] * 28313
    for fragment, bounds in [(1, (411, 15302)), (2, (14079, 28023))]:
        unit = f"{ANALYSIS}__fragment{fragment}"
        expected = audit["fragments"][f"KEL_Fragment_{fragment}"]
        summary = next(row for row in summaries if row["unit"] == unit)
        assert int(float(summary["reads"])) == expected["read_count"]
        qc = {r["metric"]: r["value"] for r in table(RESULTS / f"qc/alignment/{unit}/{unit}.alignment_qc.tsv")}
        assert int(qc["primary_records"]) == expected["read_count"]
        observed = {"unit": unit, "input_reads": expected["read_count"], "primary_mapped": int(qc["primary_mapped"])}
        for suffix in ("sorted.bam", "variant.bam", "phasing.bam"):
            path = RESULTS / f"mapping/amplicons/{unit}/{unit}.{suffix}"
            pysam.quickcheck(str(path))
            with pysam.AlignmentFile(path) as bam:
                assert bam.check_index()
                assert bam.references == ("NG_007492.3",)
                assert bam.lengths == (28313,)
                observed[suffix] = {"mapped_alignment_records": bam.mapped, "unmapped_alignment_records": bam.unmapped}
                if suffix in ("variant.bam", "phasing.bam"):
                    checked = 0
                    for read in bam.fetch(until_eof=True):
                        assert not (read.is_unmapped or read.is_secondary or read.is_supplementary)
                        assert read.mapping_quality >= 30
                        if suffix == "phasing.bam":
                            assert read.query_length >= 8000
                        checked += 1
                    assert checked == bam.mapped
                    observed[suffix]["filter_constraints_verified_reads"] = checked
        filter_log = RESULTS / f"logs/mapping/{unit}.phasing_filter.log"
        filter_stats = dict(line.split("\t", 1) for line in filter_log.read_text().splitlines())
        assert observed["phasing.bam"]["mapped_alignment_records"] == int(filter_stats["kept_records"])
        assert observed["variant.bam"]["mapped_alignment_records"] == int(qc["mapq_ge_30"])
        depth = {}
        for line in (RESULTS / f"mapping/amplicons/{unit}/{unit}.variant.depth.tsv").read_text().splitlines():
            contig, position, value = line.split("\t")
            assert contig == "NG_007492.3"
            position, value = int(position), int(value)
            depth[position] = value
            aggregate_depth[position - 1] += value
        target_depth = [depth.get(p, 0) for p in range(bounds[0], bounds[1] + 1)]
        observed["target_depth"] = depth_summary(target_depth)
        total_depth = sum(depth.values())
        outside_depth = sum(value for position, value in depth.items() if not bounds[0] <= position <= bounds[1])
        observed["fraction_depth_outside_designated_amplicon"] = outside_depth / total_depth if total_depth else None
        report["amplicons"].append(observed)
    report["combined_variant_target_depth"] = depth_summary(aggregate_depth[410:28023])
    for suffix in ("merged.bam", "variant.bam", "phasing.bam"):
        path = RESULTS / f"mapping/genes/{ANALYSIS}/{ANALYSIS}.{suffix}"
        pysam.quickcheck(str(path))
        with pysam.AlignmentFile(path) as bam:
            assert bam.check_index()
            report[suffix] = {"mapped_alignment_records": bam.mapped, "unmapped_alignment_records": bam.unmapped}
            unit_suffix = "sorted.bam" if suffix == "merged.bam" else suffix
            for metric in ("mapped_alignment_records", "unmapped_alignment_records"):
                assert report[suffix][metric] == sum(row[unit_suffix][metric] for row in report["amplicons"])
    report["status"] = "PASS"
    (SETUP / "kel_core_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k: v for k, v in report.items() if k != "sources"}, indent=2))


def depth_summary(values):
    return {
        "positions": len(values),
        "min": min(values),
        "mean": sum(values) / len(values),
        "zero_depth": sum(v == 0 for v in values),
        "below_20": sum(v < 20 for v in values),
        "below_50": sum(v < 50 for v in values),
    }


if __name__ == "__main__":
    main()
