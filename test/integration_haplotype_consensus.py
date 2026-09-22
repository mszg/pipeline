#!/usr/bin/env python3
"""
Exercise support, coverage, masking and consensus against synthetic known truth.

Run with Python, pysam, samtools, bcftools and tabix available in the environment:
    python test/integration_haplotype_consensus.py
    python test/integration_haplotype_consensus.py --output /path/to/empty/directory

Without --output, all generated fixtures are removed after the test.
"""
import argparse
import csv
import json
import random
import runpy
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from types import SimpleNamespace as NS

import pysam

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
from check_phase_sets import check_phase_sets


def run_script(name, **kwargs):
    runpy.run_path(str(SCRIPTS / name), init_globals={"snakemake": NS(**{k: NS(**v) if isinstance(v, dict) else v for k, v in kwargs.items()})})


def shell(args, **kwargs):
    return subprocess.run(args, check=True, text=True, **kwargs)


def run_integration(out):
    rng = random.Random(20260921)
    ref = list("".join(rng.choice("ACGT") for _ in range(500)))
    for start, sequence in [(50, "A"), (100, "CG"), (200, "ACG"), (300, "ACG"), (400, "A"), (450, "T")]:
        ref[start:start + len(sequence)] = sequence
    ref = "".join(ref)
    fasta = out / "reference.fasta"
    fasta.write_text(">SYNTHETIC_REF\n" + ref + "\n")
    pysam.faidx(str(fasta))
    bed = out / "targets.bed"
    bed.write_text("SYNTHETIC_REF\t10\t490\n")
    variants = [(51, "A", "G", "0|1"), (101, "C", "CT", "0|1"),
                (201, "AC", "A", "0|1"), (301, "AC", "A", "1/1"),
                (401, "A", "C,G", "1/2"), (451, "T", "C", "0/1")]
    header = '##fileformat=VCFv4.2\n##contig=<ID=SYNTHETIC_REF,length=500>\n##FILTER=<ID=PASS,Description="All filters passed">\n'
    header += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
    header += '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Depth">\n##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allele depths">\n'
    header += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSYNTHETIC__INDEL\n"
    full, phased = out / "full.vcf.gz", out / "phased.vcf.gz"
    for path, include_multi in [(full, True), (phased, False)]:
        plain = path.with_suffix("")
        rows = []
        for pos, ra, aa, gt in variants:
            if "," in aa and not include_multi:
                continue
            rows.append(f"SYNTHETIC_REF\t{pos}\t.\t{ra}\t{aa}\t50\tPASS\t.\tGT:PS:DP:AD\t{gt}:{51 if '|' in gt else '.'}:120:{'0,60,60' if ',' in aa else '60,60'}\n")
        plain.write_text(header + "".join(rows))
        pysam.tabix_compress(str(plain), str(path), force=True)
        pysam.tabix_index(str(path), preset="vcf", force=True)
    check_phase_sets(phased, bed, "SYNTHETIC__INDEL")
    bams = {}
    for hp in (1, 2):
        path = out / f"HP{hp}.bam"
        bams[hp] = path
        with pysam.AlignmentFile(path, "wb", header={"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "SYNTHETIC_REF", "LN": 500}]}) as bam:
            for index in range(60):
                stop = 500 if hp == 1 else 480
                pieces, cigar, cursor = [], [], 0
                for pos, ra, aa, gt in variants:
                    start = pos - 1
                    if start > cursor:
                        pieces.append(ref[cursor:start]); cigar.append((0, start - cursor))
                    observed = ra
                    if pos == 301:  # old parser falsely calls this exact ALT
                        observed = "G"
                    elif pos == 401:
                        observed = "C" if hp == 1 else "G"
                    elif pos == 451:
                        observed = "C" if index < 30 else "T"
                    elif hp == 2:
                        observed = aa
                    pieces.append(observed)
                    shared = min(len(ra), len(observed))
                    cigar.append((0, shared))
                    if len(ra) > len(observed):
                        cigar.append((2, len(ra) - len(observed)))
                    elif len(observed) > len(ra):
                        cigar.append((1, len(observed) - len(ra)))
                    cursor = start + len(ra)
                pieces.append(ref[cursor:stop]); cigar.append((0, stop - cursor))
                record = pysam.AlignedSegment()
                record.query_name = f"synthetic_hp{hp}_{index}"
                record.reference_id = 0
                record.reference_start = 0
                record.mapping_quality = 60
                record.flag = 0
                record.query_sequence = "".join(pieces)
                record.cigartuples = cigar
                record.query_qualities = pysam.qualitystring_to_array("I" * len(record.query_sequence))
                record.set_tag("HP", hp); record.set_tag("PS", 51)
                bam.write(record)
        pysam.index(str(path))
    outputs = {name: str(out / filename) for name, filename in {
        "support": "support.tsv", "uncertain": "uncertain.tsv", "uncertain_bed": "uncertain.bed",
        "accepted_bed": "accepted.bed", "consensus_vcf": "consensus_ready.vcf.gz", "consensus_tbi": "consensus_ready.vcf.gz.tbi"}.items()}
    params = dict(min_support_depth=20, min_het_alt_fraction=0.3, min_het_delta=0.25,
                  min_hom_alt_fraction=0.8, max_other_fraction=0.25, require_pass=True, mpileup_max_depth=100000)
    run_script("haplotype_variant_support.py", input=dict(full_vcf=str(full), phased_vcf=str(phased), hp1_bam=str(bams[1]), hp2_bam=str(bams[2]), ref=str(fasta)), output=outputs, params=params)
    with open(outputs["support"]) as handle:
        support = {int(row["POS"]): row for row in csv.DictReader(handle, delimiter="\t")}
    assert {p for p, row in support.items() if row["STATUS"] == "ACCEPT"} == {51, 101, 201}
    assert support[301]["HP1_ALT"] == support[301]["HP2_ALT"] == "0"
    assert support[301]["HP1_OTHER"] == support[301]["HP2_OTHER"] == "60"
    run_script("haplotype_coverage_qc.py", input=dict(hp1_bam=str(bams[1]), hp2_bam=str(bams[2]), bed=str(bed)), output=dict(summary=str(out / "coverage.tsv"), hp1_low=str(out / "HP1.low.bed"), hp2_low=str(out / "HP2.low.bed")), params=dict(min_depth=50))
    result = {"scope": "Synthetic exact-allele, uncertainty and consensus integration; not assay accuracy", "support_statuses": {str(p): r["STATUS"] for p, r in support.items()}, "consensus": {}}
    for hp in (1, 2):
        mask = out / f"HP{hp}.mask.bed"
        run_script("make_consensus_mask.py", input=dict(fai=str(fasta) + ".fai", target_bed=str(bed), uncertain_bed=outputs["uncertain_bed"], low_depth_bed=str(out / f"HP{hp}.low.bed"), support_tsv=outputs["support"]), output=dict(bed=str(mask)), params=dict(haplotype=hp, callable_min_depth=50))
        consensus = out / f"HP{hp}.fasta"
        with consensus.open("w") as handle:
            shell(["bcftools", "consensus", "-f", str(fasta), "-H", str(hp), "-m", str(mask), outputs["consensus_vcf"]], stdout=handle)
        expected = list(ref)
        for i in list(range(10)) + list(range(490 if hp == 1 else 480, 500)) + [300, 301, 400, 450]:
            expected[i] = "N"
        if hp == 2:
            expected[50] = "G"
            expected[200:202] = list("A")
            expected[100:101] = list("CT")
        actual = "".join(consensus.read_text().splitlines()[1:])
        result["consensus"][f"HP{hp}"] = {"length": len(actual), "expected_length": len(expected), "exact_sequence_match": actual == "".join(expected)}
    result["blocked_deletion_anchor"] = blocked_deletion_anchor(out / "blocked_deletion_anchor")
    result["status"] = "PASS" if all(x["exact_sequence_match"] for x in result["consensus"].values()) and result["blocked_deletion_anchor"]["exact_sequence_match"] else "FAIL"
    (out / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    assert result["status"] == "PASS", "Consensus differs from known synthetic haplotype/mask truth"
    return result


def blocked_deletion_anchor(out):
    """A masked anchor must prevent rescue of a deletion bcftools will skip."""
    out.mkdir()
    reference = "G" * 10 + "AC" + "G" * 8
    fasta = out / "reference.fasta"
    fasta.write_text(">BLOCKED_REF\n" + reference + "\n")
    pysam.faidx(str(fasta))
    target = out / "target.bed"
    target.write_text("BLOCKED_REF\t0\t20\n")
    uncertain = out / "uncertain.bed"
    uncertain.write_text("BLOCKED_REF\t10\t11\n")
    low = out / "low.bed"
    low.write_text("BLOCKED_REF\t11\t12\n")
    support = out / "support.tsv"
    support_row = dict(CHROM="BLOCKED_REF", POS="11", REF="AC", ALT="A",
                       TYPE="DEL", STATUS="ACCEPT", PHASED_GT="0|1",
                       HP1_ALT="0", HP2_ALT="60")
    with support.open("w") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(support_row), delimiter="\t")
        writer.writeheader()
        writer.writerow(support_row)
    plain = out / "accepted.vcf"
    plain.write_text(
        '##fileformat=VCFv4.2\n##contig=<ID=BLOCKED_REF,length=20>\n'
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tTEST\n'
        'BLOCKED_REF\t11\t.\tAC\tA\t60\tPASS\t.\tGT\t0|1\n'
    )
    vcf = out / "accepted.vcf.gz"
    pysam.tabix_compress(str(plain), str(vcf), force=True)
    pysam.tabix_index(str(vcf), preset="vcf", force=True)
    mask = out / "mask.bed"
    run_script(
        "make_consensus_mask.py",
        input=dict(fai=str(fasta) + ".fai", target_bed=str(target),
                   uncertain_bed=str(uncertain), low_depth_bed=str(low),
                   support_tsv=str(support)),
        output=dict(bed=str(mask)),
        params=dict(haplotype=2, callable_min_depth=50),
    )
    masked_positions = {
        position
        for line in mask.read_text().splitlines()
        for _, start, end in [line.split("\t")]
        for position in range(int(start), int(end))
    }
    assert masked_positions == {10, 11}, "Blocked deletion leaked rescued reference positions"
    result = shell(["bcftools", "consensus", "-f", str(fasta), "-H", "2",
                    "-m", str(mask), str(vcf)], capture_output=True)
    (out / "consensus.fasta").write_text(result.stdout)
    (out / "consensus.log").write_text(result.stderr)
    actual = "".join(result.stdout.splitlines()[1:])
    expected = "G" * 10 + "NN" + "G" * 8
    assert actual == expected, "Masked deletion anchor exposed an unsupported reference allele"
    return {"exact_sequence_match": actual == expected, "masked_positions_0based": sorted(masked_positions),
            "length": len(actual), "expected_length": len(expected)}


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output", type=Path, help="Preserve fixtures in this new or empty directory")
    args = parser.parse_args()
    missing = [tool for tool in ("samtools", "bcftools", "tabix") if shutil.which(tool) is None]
    if missing:
        parser.error("Required tools are missing from PATH: " + ", ".join(missing))
    if args.output is not None:
        out = args.output.expanduser().resolve()
        if out.exists() and (not out.is_dir() or any(out.iterdir())):
            parser.error("--output must be a new or empty directory")
        out.mkdir(parents=True, exist_ok=True)
        result = run_integration(out)
        result["output_directory"] = str(out)
    else:
        with tempfile.TemporaryDirectory(prefix="haplotype_consensus_integration_") as directory:
            result = run_integration(Path(directory))
        result["output_directory"] = None
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
