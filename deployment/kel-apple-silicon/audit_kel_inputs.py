#!/usr/bin/env python3
"""Read-only KEL FASTQ audit; no read identifiers or sequence data are emitted.

Uses only the Python standard library. Full gzip consumption validates CRC/trailers.
SHA-256 covers the original compressed bytes. Writes results only to --output.
"""
import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys

METADATA_KEYS = (
    "barcode", "barcode_alias", "basecall_model_version_id", "flow_cell_id",
    "runid", "protocol_group_id", "sample_id", "basecall_gpu",
)


def length_stats(hist):
    n = sum(hist.values())
    bases = sum(length * count for length, count in hist.items())
    if not n:
        return {"read_count": 0, "total_bases": 0}
    result = {"read_count": n, "total_bases": bases,
              "min_length": min(hist), "max_length": max(hist),
              "mean_length": bases / n}
    thresholds = {name: max(1, math.ceil(q * n)) for name, q in
                  (("p10_length", .1), ("median_length", .5),
                   ("p90_length", .9), ("p95_length", .95), ("p99_length", .99))}
    cumulative = 0
    for length, count in sorted(hist.items()):
        cumulative += count
        for name, threshold in thresholds.items():
            if name not in result and cumulative >= threshold:
                result[name] = length
    cumulative = 0
    for length, count in sorted(hist.items(), reverse=True):
        cumulative += length * count
        if cumulative >= bases / 2:
            result["n50_length"] = length
            break
    result["length_bins"] = {
        "lt_1000": sum(c for l, c in hist.items() if l < 1000),
        "1000_to_4999": sum(c for l, c in hist.items() if 1000 <= l < 5000),
        "5000_to_9999": sum(c for l, c in hist.items() if 5000 <= l < 10000),
        "10000_to_19999": sum(c for l, c in hist.items() if 10000 <= l < 20000),
        "ge_20000": sum(c for l, c in hist.items() if l >= 20000),
    }
    return result


def quality_stats(hist):
    n = sum(hist.values())
    if not n:
        return {}
    mean_error = sum(c * 10 ** (-(q - 33) / 10) for q, c in hist.items()) / n
    return {
        "quality_encoding_assumption": "Phred+33",
        "min_base_q": min(hist) - 33,
        "max_base_q": max(hist) - 33,
        "arithmetic_mean_base_q": sum((q - 33) * c for q, c in hist.items()) / n,
        "q_from_mean_base_error_probability": -10 * math.log10(mean_error),
        "base_fraction_q_ge_10": sum(c for q, c in hist.items() if q >= 43) / n,
        "base_fraction_q_ge_20": sum(c for q, c in hist.items() if q >= 53) / n,
        "base_fraction_q_ge_30": sum(c for q, c in hist.items() if q >= 63) / n,
    }


def scan(path, source):
    before = path.stat()
    lengths, qualities = collections.Counter(), collections.Counter()
    metadata = {key: collections.Counter() for key in METADATA_KEYS}
    result = {"relative_path": str(path.relative_to(source)),
              "compressed_bytes": before.st_size,
              "mtime_ns": before.st_mtime_ns,
              "gzip_integrity": "pending", "fastq_integrity": "pending"}
    times = []
    first_time = last_time = None
    missing_metadata = collections.Counter()
    invalid_nucleotides = collections.Counter()
    noncanonical_bases = 0
    record = 0
    try:
        with gzip.open(path, "rb") as stream:
            while True:
                header = stream.readline()
                if not header:
                    break
                record += 1
                sequence = stream.readline().rstrip(b"\r\n")
                plus = stream.readline()
                quality = stream.readline().rstrip(b"\r\n")
                if not header.startswith(b"@") or not plus.startswith(b"+"):
                    raise ValueError(f"Invalid FASTQ record structure at record {record}")
                if not sequence or len(sequence) != len(quality):
                    raise ValueError(f"Empty sequence or length mismatch at record {record}")
                if min(quality) < 33 or max(quality) > 126:
                    raise ValueError(f"Invalid Phred+33 quality byte at record {record}")
                invalid = sequence.translate(None, b"ACGTNacgtn")
                if invalid:
                    invalid_nucleotides.update(invalid)
                    noncanonical_bases += len(invalid)
                lengths[len(sequence)] += 1
                qualities.update(quality)
                fields = {}
                for item in header.rstrip(b"\r\n").split()[1:]:
                    if b"=" in item:
                        key, value = item.split(b"=", 1)
                        fields[key.decode("ascii")] = value.decode("utf-8")
                for key in METADATA_KEYS:
                    if key in fields:
                        metadata[key][fields[key]] += 1
                    else:
                        missing_metadata[key] += 1
                value = fields.get("start_time")
                if value:
                    first_time = value if first_time is None else min(first_time, value)
                    last_time = value if last_time is None else max(last_time, value)
        result["gzip_integrity"] = "pass"
        result["fastq_integrity"] = "pass" if not noncanonical_bases else "non_ACGTN_bases"
    except Exception as exc:
        # Exceptions do not contain raw read data.
        result["error"] = f"{type(exc).__name__}: {exc}"
        result["gzip_integrity"] = "failed_or_incomplete"
        result["fastq_integrity"] = "failed_or_incomplete"
    with path.open("rb") as stream:
        digest = hashlib.file_digest(stream, "sha256").hexdigest()
    after = path.stat()
    result.update(length_stats(lengths))
    result.update(quality_stats(qualities))
    result.update({"sha256": digest,
                   "source_size_and_mtime_unchanged":
                       (before.st_size, before.st_mtime_ns) == (after.st_size, after.st_mtime_ns),
                   "header_metadata_value_counts": metadata,
                   "missing_header_metadata_record_counts": missing_metadata,
                   "non_ACGTN_base_count": noncanonical_bases,
                   "first_read_start_time": first_time,
                   "last_read_start_time": last_time})
    return result, lengths, qualities, metadata


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    source = args.source.resolve()
    output = args.output.resolve()
    if source == output or source in output.parents:
        parser.error("Output must be outside the original raw data tree")
    output.mkdir(parents=True, exist_ok=True)
    paths = sorted(source.rglob("*.fastq.gz"), key=lambda p: (p.parent.name, int(re.search(r"_(\d+)\.fastq\.gz$", p.name).group(1))))
    report = {"source_directory": str(source),
              "audit_started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python_version": sys.version,
              "method": "Read-only full gzip and 4-line FASTQ validation, SHA-256 of compressed originals; aggregate Phred+33 qualities; nearest-rank length quantiles; no sample identity inference",
              "files": [], "fragments": {}, "status": "in_progress"}
    aggregate = {}
    for path in paths:
        result, lengths, qualities, metadata = scan(path, source)
        report["files"].append(result)
        fragment = path.parent.name
        if fragment not in aggregate:
            aggregate[fragment] = {"lengths": collections.Counter(), "qualities": collections.Counter(),
                                   "metadata": {key: collections.Counter() for key in METADATA_KEYS},
                                   "indices": [], "compressed_bytes": 0}
        agg = aggregate[fragment]
        agg["lengths"].update(lengths)
        agg["qualities"].update(qualities)
        for key in METADATA_KEYS:
            agg["metadata"][key].update(metadata[key])
        agg["indices"].append(int(re.search(r"_(\d+)\.fastq\.gz$", path.name).group(1)))
        agg["compressed_bytes"] += result["compressed_bytes"]
        (output / "kel_input_audit.partial.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps({"completed_files": len(report["files"]), "total_files": len(paths),
                          "file": result["relative_path"], "reads": result.get("read_count", 0),
                          "gzip_integrity": result["gzip_integrity"], "fastq_integrity": result["fastq_integrity"]}), flush=True)
    for fragment, agg in aggregate.items():
        indices = sorted(agg["indices"])
        report["fragments"][fragment] = {
            "file_count": len(indices), "compressed_bytes": agg["compressed_bytes"],
            "observed_chunk_indices": indices,
            "absent_indices_within_observed_range": sorted(set(range(min(indices), max(indices) + 1)) - set(indices)),
            "chunk_gap_interpretation": "A filename gap does not prove missing sequencing data; verify against the instrument export/sample sheet",
            "header_metadata_value_counts": agg["metadata"],
            **length_stats(agg["lengths"]), **quality_stats(agg["qualities"]),
        }
        with (output / f"{fragment}.length_histogram.tsv").open("w") as stream:
            stream.write("read_length_bp\tread_count\n")
            for length, count in sorted(agg["lengths"].items()):
                stream.write(f"{length}\t{count}\n")
        with (output / f"{fragment}.base_quality_histogram.tsv").open("w") as stream:
            stream.write("phred_q\tbase_count\n")
            for quality, count in sorted(agg["qualities"].items()):
                stream.write(f"{quality - 33}\t{count}\n")
    report["audit_finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
    report["status"] = "pass" if all(x["gzip_integrity"] == "pass" and x["fastq_integrity"] == "pass" and x["source_size_and_mtime_unchanged"] for x in report["files"]) else "review_required"
    (output / "kel_input_audit.json").write_text(json.dumps(report, indent=2) + "\n")
    columns = ["relative_path", "compressed_bytes", "sha256", "read_count", "total_bases",
               "min_length", "median_length", "p90_length", "max_length", "n50_length",
               "arithmetic_mean_base_q", "q_from_mean_base_error_probability",
               "gzip_integrity", "fastq_integrity", "source_size_and_mtime_unchanged"]
    with (output / "kel_input_inventory.tsv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(report["files"])
    print(json.dumps({"status": report["status"], "fragments": report["fragments"]}), flush=True)
    return 0 if report["status"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
