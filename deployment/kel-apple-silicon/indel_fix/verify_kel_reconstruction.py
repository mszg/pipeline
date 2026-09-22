#!/usr/bin/env python3
"""Independent KEL output accounting and sequence reconstruction checks."""
import collections
import csv
import hashlib
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parents[2]
RESULTS = ROOT / "workflow/results"
A = "KEL__KEL"
CHROM = "NG_007492.3"
QC = RESULTS / "qc/haplotypes" / A
HAPS = RESULTS / "haplotypes" / A
CONS = RESULTS / "consensus" / A


def rows(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def positions(path):
    result = set()
    for line in path.read_text().splitlines():
        chrom, start, end, *_ = line.split("\t")
        assert chrom == CHROM
        result.update(range(int(start), int(end)))
    return result


def key(row):
    return row["CHROM"], int(row["POS"]), row["REF"], row["ALT"]


def main():
    report = {"scope": "Read partition, support accounting, mask and complete consensus consistency; no external biological truth"}
    names = {1: set(), 2: set(), 0: set()}
    with pysam.AlignmentFile(HAPS / f"{A}.haplotagged.bam", "rb") as bam:
        assert bam.check_index()
        for read in bam:
            hp = read.get_tag("HP") if read.has_tag("HP") else 0
            assert hp in names and read.query_name not in names[hp]
            if hp:
                assert read.get_tag("PS") == 516
            names[hp].add(read.query_name)
    assert not (names[1] & names[2] or names[1] & names[0] or names[2] & names[0])
    with pysam.AlignmentFile(RESULTS / "mapping/genes" / A / f"{A}.phasing.bam", "rb") as bam:
        original_names = {r.query_name for r in bam}
    assert set.union(*names.values()) == original_names and len(original_names) == 85890
    for hp in (1, 2):
        with pysam.AlignmentFile(HAPS / f"{A}.HP{hp}.bam", "rb") as bam:
            assert bam.check_index()
            observed = set()
            for read in bam:
                assert read.get_tag("HP") == hp and read.get_tag("PS") == 516
                assert read.query_name not in observed
                observed.add(read.query_name)
            assert observed == names[hp]
    report["read_partition"] = {f"HP{hp}" if hp else "unassigned": len(value) for hp, value in names.items()}
    phase = rows(QC / f"{A}.phase_set_validation.tsv")
    assert len(phase) == 1 and phase[0]["phase_set"] == "516" and phase[0]["status"] == "PASS"
    support = rows(QC / f"{A}.variant_support.tsv")
    accepted = [r for r in support if r["STATUS"] == "ACCEPT"]
    unresolved = [r for r in support if r["STATUS"] == "UNRESOLVED"]
    assert len(support) == 45 and len(accepted) + len(unresolved) == len(support)
    assert rows(QC / f"{A}.uncertain_variants.tsv") == unresolved
    with pysam.VariantFile(RESULTS / "variants" / A / f"{A}.norm.vcf.gz") as vcf:
        assert {key(r) for r in support} == {(r.contig, r.pos, r.ref, ",".join(r.alts)) for r in vcf.fetch()}
    with pysam.VariantFile(RESULTS / "variants" / A / f"{A}.consensus_ready.vcf.gz") as vcf:
        assert list(vcf.header.samples) == [A]
        records = list(vcf.fetch())
        assert {key(r) for r in accepted} == {(r.contig, r.pos, r.ref, ",".join(r.alts)) for r in records}
        by_key = {key(r): r for r in accepted}
        for rec in records:
            row = by_key[(rec.contig, rec.pos, rec.ref, rec.alts[0])]
            sample = rec.samples[A]
            assert tuple(map(int, row["PHASED_GT"].replace("|", "/").split("/"))) == sample["GT"]
            assert sample["GT"][0] == sample["GT"][1] or sample.phased
    uncertain = {p for r in unresolved for p in range(int(r["POS"]) - 1, int(r["POS"]) - 1 + len(r["REF"]))}
    assert uncertain == positions(QC / f"{A}.uncertain_regions.bed")
    report["status_counts"] = dict(collections.Counter(r["STATUS"] for r in support))
    report["accepted_type_counts"] = dict(collections.Counter(r["TYPE"] for r in accepted))
    report["unresolved"] = [{c: r[c] for c in ("POS", "REF", "ALT", "TYPE", "REASON")} for r in unresolved]
    # Independently verify the indel denominator and sample classifications using
    # aligned-pair query coordinates (not the production CIGAR-walking extractor).
    import sys
    sys.path.insert(0, str(ROOT / "workflow/workflow/scripts"))
    from allele_support import indel_window, classify_indel
    reference = "".join((RESULTS / "reference" / A / "reference.fasta").read_text().splitlines()[1:])
    sample_checks = 0
    for hp in (1, 2):
        with pysam.AlignmentFile(HAPS / f"{A}.HP{hp}.bam", "rb") as bam:
            for row in support:
                if row["TYPE"] not in {"INS", "DEL"}:
                    continue
                start = int(row["POS"]) - 1
                window = indel_window(reference, start, row["REF"], row["ALT"])
                count = 0
                for read in bam.fetch(CHROM, start, start + 1):
                    if read.flag & (4 | 256 | 512 | 1024 | 2048):
                        continue
                    count += 1
                    if count <= 20:
                        pairs = dict((r, q) for q, r in read.get_aligned_pairs() if r is not None)
                        left, right = pairs.get(window.start), pairs.get(window.end - 1)
                        observed = read.query_sequence[left:right + 1] if left is not None and right is not None else None
                        expected = "REF" if observed == window.ref else "ALT" if observed == window.alt else "OTHER"
                        # KEL is DNA alignment: no reference skips or internal clips.
                        assert all(op not in (3, 6) for op, _ in read.cigartuples)
                        assert classify_indel(read, window) == expected
                        sample_checks += 1
                assert count == sum(int(row[f"HP{hp}_{c}"]) for c in ("REF", "ALT", "OTHER"))
    report["indel_independent_sample_checks"] = sample_checks
    target = positions(RESULTS / "targets" / A / f"{A}.bed")
    assert target == set(range(410, 28023))
    outside = set(range(len(reference))) - target
    report["consensus"] = {}
    for hp in (1, 2):
        low = positions(QC / f"{A}.HP{hp}.low_depth.bed")
        deletion_bases = set()
        for row in accepted:
            allele = row["PHASED_GT"].replace("|", "/").split("/")[hp - 1]
            if row["TYPE"] == "DEL" and allele == "1" and int(row[f"HP{hp}_ALT"]) >= 50:
                start = int(row["POS"]) - 1
                span = set(range(start, start + len(row["REF"])))
                if not (span & (outside | uncertain)) and start not in low:
                    deletion_bases.update(span - {start})
        mask = positions(CONS / f"{A}.HP{hp}.mask.bed")
        assert mask == outside | uncertain | (low - deletion_bases)
        expected = list(reference)
        for p in mask:
            expected[p] = "N"
        applied, skipped = [], []
        for row in sorted(accepted, key=lambda r: int(r["POS"]), reverse=True):
            if row["PHASED_GT"].replace("|", "/").split("/")[hp - 1] != "1":
                continue
            start, ref, alt = int(row["POS"]) - 1, row["REF"], row["ALT"]
            assert reference[start:start + len(ref)] == ref
            if set(range(start, start + len(ref))) & mask:
                skipped.append(int(row["POS"]))
                continue
            expected[start:start + len(ref)] = list(alt)
            applied.append(int(row["POS"]))
        actual = "".join((CONS / f"{A}.haplotype{hp}.fasta").read_text().splitlines()[1:])
        assert actual == "".join(expected), f"HP{hp} sequence mismatch"
        report["consensus"][f"HP{hp}"] = {"length": len(actual), "N_bases": actual.count("N"), "target_masked_reference_bases": len(mask & target), "target_unmasked_reference_bases": len(target - mask), "low_depth_reference_bases": len(low), "rescued_deletion_reference_bases": len((low & deletion_bases) - uncertain - outside), "applied_alt_positions": sorted(applied), "accepted_alt_positions_masked": sorted(skipped), "exact_independent_sequence_match": True, "sha256": hashlib.sha256(actual.encode()).hexdigest()}
    audit = json.loads((ROOT / "setup/input_audit/kel_input_audit.json").read_text())
    for entry in audit["files"]:
        path = Path(audit["source_directory"]) / entry["relative_path"]
        with path.open("rb") as handle:
            digest = hashlib.file_digest(handle, "sha256").hexdigest()
        assert digest == entry["sha256"], path
    report["original_files_unchanged"] = len(audit["files"])
    report["status"] = "PASS"
    (Path(__file__).resolve().parent / "kel_reconstruction_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
