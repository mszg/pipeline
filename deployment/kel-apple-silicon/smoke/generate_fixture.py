#!/usr/bin/env python3
"""Generate deterministic synthetic data to exercise the core workflow tools."""

import csv
import gzip
import hashlib
import json
import random
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SEED = 20260921
CONTIG = "SYNTHETIC_REF"
BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGT", "TGCA")
AMPLICON_BOUNDS = {"amplicon1": (1001, 14000), "amplicon2": (10001, 23000)}


def main():
    rng = random.Random(SEED)
    reference = "".join(rng.choice(BASES) for _ in range(24000))
    reference_path = ROOT / "reference.fasta"
    with reference_path.open("w") as handle:
        handle.write(f">{CONTIG}\n")
        for i in range(0, len(reference), 80):
            handle.write(reference[i:i + 80] + "\n")

    generated = [reference_path]
    expected_reads = []
    sample_rows = []
    for amplicon, (start, end) in AMPLICON_BOUNDS.items():
        input_dir = ROOT / "fastqs" / amplicon
        input_dir.mkdir(parents=True, exist_ok=True)
        records = []
        for i in range(24):
            # Four full-span reads ensure deterministic target-edge coverage.
            left_trim = 0 if i < 4 else rng.randrange(201)
            right_trim = 0 if i < 4 else rng.randrange(201)
            read_start, read_end = start + left_trim, end - right_trim
            sequence = list(reference[read_start - 1:read_end])
            # Substitutions leave a simple known alignment span; protect ends.
            substitutions = round(len(sequence) * (0.004 + (i % 5) * 0.001))
            for position in rng.sample(range(50, len(sequence) - 50), substitutions):
                sequence[position] = rng.choice(BASES.replace(sequence[position], ""))
            sequence = "".join(sequence)
            qualities = "".join(chr(33 + 20 + (i % 9) + rng.randrange(4)) for _ in sequence)
            strand = "+" if i % 2 == 0 else "-"
            if strand == "-":
                sequence = sequence.translate(COMPLEMENT)[::-1]
                qualities = qualities[::-1]
            read_id = f"SYNTHETIC_{amplicon}_{i + 1:03d}"
            records.append(f"@{read_id}\n{sequence}\n+\n{qualities}\n")
            expected_reads.append({
                "read_id": read_id,
                "amplicon": amplicon,
                "contig": CONTIG,
                "start_1based": read_start,
                "end_1based": read_end,
                "strand": strand,
                "length": len(sequence),
                "substitutions": substitutions,
            })

        for chunk in range(3):
            suffix = ".fastq.gz" if chunk == 1 else ".fastq"
            path = input_dir / f"chunk_{chunk + 1:02d}{suffix}"
            content = "".join(records[chunk * 8:(chunk + 1) * 8]).encode("ascii")
            if suffix.endswith(".gz"):
                with path.open("wb") as raw:
                    with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as handle:
                        handle.write(content)
            else:
                path.write_bytes(content)
            generated.append(path)
        sample_rows.append([
            "SYNTHETIC", "TEST", amplicon, str(input_dir), str(reference_path),
            f"{CONTIG}:{start}-{end}", "8000", "3",
        ])

    with (ROOT / "samples.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "sample", "gene", "amplicon", "fastq_input", "reference",
            "target_region", "phasing_min_length", "core_depth_threshold",
        ])
        writer.writerows(sample_rows)

    with (ROOT / "expected_reads.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(expected_reads[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(expected_reads)

    expected = {
        "purpose": "Synthetic core workflow execution smoke test; not biological validation",
        "seed": SEED,
        "reference_contig": CONTIG,
        "reference_length": len(reference),
        "analysis": "SYNTHETIC__TEST",
        "amplicons": 2,
        "fastq_chunks": 6,
        "reads_per_chunk": 8,
        "reads_per_amplicon": 24,
        "total_reads": len(expected_reads),
        "total_bases": sum(row["length"] for row in expected_reads),
        "min_read_length": min(row["length"] for row in expected_reads),
        "max_read_length": max(row["length"] for row in expected_reads),
        "expected_primary_mapped_reads": len(expected_reads),
        "expected_variant_reads": len(expected_reads),
        "expected_phasing_reads": len(expected_reads),
        "expected_merged_target_bed": f"{CONTIG}\t1000\t23000",
        "target_union_bases": 22000,
    }
    (ROOT / "expected.json").write_text(json.dumps(expected, indent=2) + "\n")
    (ROOT / "expected.targets.bed").write_text(expected["expected_merged_target_bed"] + "\n")
    with (ROOT / "sha256.tsv").open("w") as handle:
        handle.write("sha256\tfile\n")
        for path in sorted(generated):
            handle.write(f"{hashlib.sha256(path.read_bytes()).hexdigest()}\t{path.relative_to(ROOT)}\n")
    print(json.dumps(expected, indent=2))


if __name__ == "__main__":
    main()
