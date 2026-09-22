#!/usr/bin/env python3
"""Withhold individual KEL variants from haplotagging, then reevaluate support.

This is a diagnostic conditional on the existing phase orientation. The original
phasing solution and genotype remain fixed; this is not independent phase truth.
Only scratch files under leave_one_out/ are written.
"""
import csv
import hashlib
import json
import os
import runpy
import shlex
import subprocess
import sys
import time
from collections import Counter
from pathlib import Path
from types import SimpleNamespace as NS

import pysam


ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "leave_one_out"
OUT.mkdir(exist_ok=True)
ENV = ROOT / "setup/envs/kel-native/bin"
os.environ["PATH"] = str(ENV) + os.pathsep + os.environ.get("PATH", "")
for variable in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[variable] = "1"
SCRIPTS = ROOT / "workflow/workflow/scripts"
sys.path.insert(0, str(SCRIPTS))
A = "KEL__KEL"
R = ROOT / "workflow/results"
FULL = R / "variants" / A / f"{A}.norm.vcf.gz"
PHASED = R / "phasing" / A / f"{A}.phased.vcf.gz"
BAM = R / "mapping/genes" / A / f"{A}.phasing.bam"
TAGGED = R / "haplotypes" / A / f"{A}.haplotagged.bam"
REF = R / "reference" / A / "reference.fasta"
SUPPORT = R / "qc/haplotypes" / A / f"{A}.variant_support.tsv"
CASES = [("NG_007492.3", 560, "CCT", "C"), ("NG_007492.3", 714, "G", "C")]
PARAMS = dict(min_support_depth=20, min_het_alt_fraction=.30, min_het_delta=.25,
              min_hom_alt_fraction=.80, max_other_fraction=.25, require_pass=True,
              mpileup_max_depth=100000)


def sha256(path):
    with path.open("rb") as handle:
        return hashlib.file_digest(handle, "sha256").hexdigest()


def read_support(path):
    with path.open() as handle:
        return {(r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]): r
                for r in csv.DictReader(handle, delimiter="\t")}


def assignments(path):
    result = {}
    phase_sets = Counter()
    with pysam.AlignmentFile(path, "rb") as bam:
        for read in bam.fetch(until_eof=True):
            key = (read.query_name, read.flag, read.reference_id, read.reference_start)
            assert key not in result, f"Duplicate alignment key: {key}"
            hp = str(read.get_tag("HP")) if read.has_tag("HP") else "unassigned"
            result[key] = hp
            if read.has_tag("PS"):
                phase_sets[str(read.get_tag("PS"))] += 1
    return result, dict(phase_sets)


def run(command, destination):
    command = [str(argument) for argument in command]
    started = time.monotonic()
    with (destination / "commands.log").open("a") as log:
        log.write("$ " + shlex.join(command) + "\n")
        log.flush()
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
        log.write(f"[completed in {time.monotonic() - started:.3f} seconds]\n")


before_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in (FULL, PHASED, SUPPORT)}
original = read_support(SUPPORT)
original_assignments, original_ps = assignments(TAGGED)
assert len(original_assignments) == 85890
results = []
for variant in CASES:
    chrom, pos, ref, alt = variant
    case = OUT / str(pos)
    case.mkdir(exist_ok=True)
    print(f"Starting leave-one-out position {pos}", flush=True)
    scratch_vcf = case / "tagging_without_site.vcf.gz"
    removed = []
    with pysam.VariantFile(PHASED) as source:
        with pysam.VariantFile(scratch_vcf, "wz", header=source.header) as dest:
            for record in source:
                key = (record.contig, record.pos, record.ref, ",".join(record.alts))
                if key == variant:
                    removed.append(str(record).rstrip())
                else:
                    dest.write(record)
    assert len(removed) == 1, removed
    (case / "withheld_record.vcf.txt").write_text(removed[0] + "\n")
    run(["bcftools", "index", "-t", "-f", scratch_vcf], case)
    with pysam.VariantFile(PHASED) as source, pysam.VariantFile(scratch_vcf) as dest:
        expected_records = [str(record) for record in source
                            if (record.contig, record.pos, record.ref, ",".join(record.alts)) != variant]
        assert expected_records == [str(record) for record in dest]
    tagged = case / "haplotagged.bam"
    run(["whatshap", "haplotag", "--reference", REF, "--output", tagged, scratch_vcf, BAM], case)
    run(["samtools", "index", "-@", "2", tagged], case)
    for hp in (1, 2):
        hp_bam = case / f"HP{hp}.bam"
        run(["samtools", "view", "-@", "2", "-b", "-d", f"HP:{hp}", tagged, "-o", hp_bam], case)
        run(["samtools", "index", "-@", "2", hp_bam], case)
    output_names = dict(support="support.tsv", uncertain="uncertain.tsv", uncertain_bed="uncertain.bed",
                        accepted_bed="accepted.bed", consensus_vcf="consensus_ready.vcf.gz",
                        consensus_tbi="consensus_ready.vcf.gz.tbi")
    context = NS(input=NS(full_vcf=str(FULL), phased_vcf=str(PHASED), ref=str(REF),
                          hp1_bam=str(case / "HP1.bam"), hp2_bam=str(case / "HP2.bam")),
                 output=NS(**{name: str(case / filename) for name, filename in output_names.items()}),
                 params=NS(**PARAMS))
    with (case / "commands.log").open("a") as log:
        log.write("runpy haplotype_variant_support.py with original full/phased VCFs, scratch HP BAMs\n")
        log.write("params=" + json.dumps(PARAMS, sort_keys=True) + "\n")
    runpy.run_path(str(SCRIPTS / "haplotype_variant_support.py"), init_globals={"snakemake": context})
    observed = read_support(case / "support.tsv")
    assert original.keys() == observed.keys()
    after_assignments, phase_sets = assignments(tagged)
    assert original_assignments.keys() == after_assignments.keys()
    transitions = Counter(f"{original_assignments[key]}->{hp}" for key, hp in after_assignments.items())
    changed = [{"variant": list(key), "before": original[key]["STATUS"], "after": observed[key]["STATUS"]}
               for key in original if original[key]["STATUS"] != observed[key]["STATUS"]]
    result = dict(withheld_variant=list(variant), original=original[variant], leave_one_out=observed[variant],
                  original_assignments=dict(Counter(original_assignments.values())),
                  leave_one_out_assignments=dict(Counter(after_assignments.values())),
                  assignment_transitions=dict(sorted(transitions.items())), original_phase_sets=original_ps,
                  leave_one_out_phase_sets=phase_sets, all_variant_status_changes=changed)
    (case / "verification.json").write_text(json.dumps(result, indent=2) + "\n")
    results.append(result)
    print(json.dumps(result, indent=2), flush=True)

after_hashes = {str(path.relative_to(ROOT)): sha256(path) for path in (FULL, PHASED, SUPPORT)}
assert before_hashes == after_hashes
summary = dict(scope="Leave one variant out of haplotagging, retain existing phase orientation and original genotype for support evaluation. This tests direct tagging circularity, not independent phase truth.",
               thresholds=PARAMS, production_vcf_and_support_hashes_unchanged=before_hashes,
               cases=results)
(OUT / "verification.json").write_text(json.dumps(summary, indent=2) + "\n")
print("Both leave-one-out diagnostics completed; production VCF/support hashes unchanged.", flush=True)
