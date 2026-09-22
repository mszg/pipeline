#!/usr/bin/env python3
"""Check KEL VCF/index/reference consistency and summarize observed phase sets."""
import collections
import json
from pathlib import Path

import pysam

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "workflow/results"
ANALYSIS = "KEL__KEL"


def inspect_vcf(path, reference):
    counts = collections.Counter()
    genotypes = collections.Counter()
    phase_sets = collections.defaultdict(list)
    keys = []
    with pysam.VariantFile(path) as vcf:
        assert list(vcf.header.samples) == [ANALYSIS]
        for record in vcf.fetch():
            assert record.contig == "NG_007492.3"
            assert 411 <= record.pos <= 28023, record.pos
            assert reference.fetch(record.contig, record.start, record.start + len(record.ref)).upper() == record.ref.upper()
            keys.append((record.contig, record.pos, record.ref, record.alts))
            counts["records"] += 1
            counts["biallelic" if len(record.alts or ()) == 1 else "multiallelic"] += 1
            counts["PASS" if "PASS" in record.filter else "non_PASS"] += 1
            sample = record.samples[ANALYSIS]
            genotype = sample.get("GT", ())
            genotypes[("|" if sample.phased else "/").join("." if x is None else str(x) for x in genotype)] += 1
            if len(genotype) == 2 and None not in genotype and genotype[0] != genotype[1]:
                counts["heterozygous"] += 1
                counts["phased_heterozygous" if sample.phased else "unphased_heterozygous"] += 1
                phase_set = sample.get("PS")
                if sample.phased and phase_set is not None:
                    phase_sets[str(phase_set)].append(record.pos)
    return keys, {
        "file": str(path.relative_to(ROOT)),
        "counts": dict(counts),
        "genotypes": dict(genotypes),
        "phase_sets": {key: {"heterozygous_records": len(values), "start": min(values), "end": max(values)} for key, values in phase_sets.items()},
    }


def main():
    base = RESULTS / f"variants/{ANALYSIS}"
    with pysam.FastaFile(str(RESULTS / f"reference/{ANALYSIS}/reference.fasta")) as reference:
        norm_keys, normalized = inspect_vcf(base / f"{ANALYSIS}.norm.vcf.gz", reference)
        ready_keys, ready = inspect_vcf(base / f"{ANALYSIS}.phasing_ready.vcf.gz", reference)
        phased_keys, phased = inspect_vcf(RESULTS / f"phasing/{ANALYSIS}/{ANALYSIS}.phased.vcf.gz", reference)
    assert ready_keys == [key for key in norm_keys if len(key[3] or ()) == 1]
    assert phased_keys == ready_keys
    report = {
        "status": "PASS",
        "scope": "Indexed VCF parsing, sample/reference/target consistency, catalogue preservation and observed phase connectivity; no external variant truth set",
        "normalized": normalized,
        "phasing_ready": ready,
        "phased": phased,
        "single_observed_phase_set": len(phased["phase_sets"]) == 1,
        "downstream_validation_record": "setup/indel_fix/kel_reconstruction_verification.json",
    }
    (ROOT / "setup/kel_calling_verification.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
