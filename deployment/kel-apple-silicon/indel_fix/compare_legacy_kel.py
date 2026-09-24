#!/usr/bin/env python3
"""Run the archived original counter on identical HP BAMs, without replacing results."""
import csv
import json
import runpy
from pathlib import Path
from types import SimpleNamespace as NS

ROOT = Path(__file__).resolve().parents[2]
A = "KEL__KEL"
R = ROOT / "workflow/results"
OUT = Path(__file__).resolve().parent / "legacy_kel"
OUT.mkdir(exist_ok=True)
outputs = {name: str(OUT / filename) for name, filename in {
    "support": "support.tsv", "uncertain": "uncertain.tsv", "uncertain_bed": "uncertain.bed",
    "accepted_bed": "accepted.bed", "consensus_vcf": "consensus_ready.vcf.gz", "consensus_tbi": "consensus_ready.vcf.gz.tbi"}.items()}
context = NS(input=NS(full_vcf=str(R / "variants" / A / f"{A}.norm.vcf.gz"),
                      phased_vcf=str(R / "phasing" / A / f"{A}.phased.vcf.gz"),
                      hp1_bam=str(R / "haplotypes" / A / f"{A}.HP1.bam"),
                      hp2_bam=str(R / "haplotypes" / A / f"{A}.HP2.bam"),
                      ref=str(R / "reference" / A / "reference.fasta")),
             output=NS(**outputs),
             params=NS(min_support_depth=20, min_het_alt_fraction=.30, min_het_delta=.25,
                       min_hom_alt_fraction=.80, max_other_fraction=.25, require_pass=True, mpileup_max_depth=100000))
runpy.run_path(str(OUT.parent / "before/haplotype_variant_support.py"), init_globals={"snakemake": context})


def read(path):
    with path.open() as handle:
        return {(r["CHROM"], r["POS"], r["REF"], r["ALT"]): r for r in csv.DictReader(handle, delimiter="\t")}


before = read(OUT / "support.tsv")
after = read(R / "qc/haplotypes" / A / f"{A}.variant_support.tsv")
assert before.keys() == after.keys()
columns = ["STATUS", "REASON"] + [f"HP{hp}_{field}" for hp in (1, 2) for field in ("REF", "ALT", "OTHER", "ALT_FRAC", "OTHER_FRAC")]
differences = [{"position": int(k[1]), "ref": k[2], "alt": k[3], "type": after[k]["TYPE"], "before": {c: before[k][c] for c in columns}, "after": {c: after[k][c] for c in columns}} for k in before if any(before[k][c] != after[k][c] for c in columns)]
result = {"scope": "Archived defective parser versus corrected parser on identical real KEL HP BAMs and unchanged thresholds", "status_changes": sum(before[k]["STATUS"] != after[k]["STATUS"] for k in before), "changed_records": len(differences), "differences": differences}
(OUT.parent / "kel_legacy_comparison.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2))
