import subprocess
import tempfile
from pathlib import Path

full_vcf = str(snakemake.input.full_vcf)
phased_vcf = str(snakemake.input.phased_vcf)
hp1_bam = str(snakemake.input.hp1_bam)
hp2_bam = str(snakemake.input.hp2_bam)
ref = str(snakemake.input.ref)

support_tsv = Path(str(snakemake.output.support))
uncertain_tsv = Path(str(snakemake.output.uncertain))
uncertain_bed = Path(str(snakemake.output.uncertain_bed))
accepted_bed = Path(str(snakemake.output.accepted_bed))
consensus_vcf = Path(str(snakemake.output.consensus_vcf))
consensus_tbi = Path(str(snakemake.output.consensus_tbi))

min_support_depth = int(snakemake.params.min_support_depth)
min_het_alt_fraction = float(snakemake.params.min_het_alt_fraction)
min_het_delta = float(snakemake.params.min_het_delta)
min_hom_alt_fraction = float(snakemake.params.min_hom_alt_fraction)
max_other_fraction = float(snakemake.params.max_other_fraction)
require_pass = bool(snakemake.params.require_pass)
mpileup_max_depth = int(snakemake.params.mpileup_max_depth)

for path in [support_tsv, uncertain_tsv, uncertain_bed, accepted_bed, consensus_vcf]:
    path.parent.mkdir(parents=True, exist_ok=True)


def run(args, **kwargs):
    return subprocess.run(args, check=True, text=True, **kwargs)


def query_lines(vcf, fmt):
    result = run(["bcftools", "query", "-f", fmt, vcf], capture_output=True)
    return [line for line in result.stdout.splitlines() if line]


def variant_type(ref_allele, alt_allele):
    if "," in alt_allele:
        return "MULTIALLELIC"
    if len(ref_allele) == 1 and len(alt_allele) == 1:
        return "SNV"
    if len(ref_allele) > len(alt_allele) and ref_allele.startswith(alt_allele):
        return "DEL"
    if len(alt_allele) > len(ref_allele) and alt_allele.startswith(ref_allele):
        return "INS"
    if len(ref_allele) == len(alt_allele):
        return "MNV"
    return "COMPLEX"


def observations(bases):
    out = []
    i = 0
    while i < len(bases):
        c = bases[i]
        if c == "^":
            i += 2
            continue
        if c == "$":
            i += 1
            continue
        if c in ".,ACGTNacgtn*#<>":
            obs = {"base": c, "indels": []}
            i += 1
            while i < len(bases) and bases[i] in "+-":
                sign = bases[i]
                i += 1
                j = i
                while j < len(bases) and bases[j].isdigit():
                    j += 1
                if j == i:
                    break
                n = int(bases[i:j])
                seq = bases[j:j + n].upper()
                obs["indels"].append((sign, n, seq))
                i = j + n
            out.append(obs)
            continue
        i += 1
    return out


def count_support(ref_allele, alt_allele, bases):
    obs = observations(bases)
    ref_n = alt_n = other = 0

    if len(ref_allele) == 1 and len(alt_allele) == 1:
        for item in obs:
            base = item["base"]
            if base in ".,":
                ref_n += 1
            elif base.upper() == alt_allele.upper():
                alt_n += 1
            elif base.upper() in "ACGTN":
                other += 1
        return ref_n, alt_n, other, True

    if len(ref_allele) > len(alt_allele) and ref_allele.startswith(alt_allele):
        deleted = ref_allele[len(alt_allele):].upper()
        target = ("-", len(deleted), deleted)
    elif len(alt_allele) > len(ref_allele) and alt_allele.startswith(ref_allele):
        inserted = alt_allele[len(ref_allele):].upper()
        target = ("+", len(inserted), inserted)
    else:
        return 0, 0, len(obs), False

    for item in obs:
        if target in item["indels"]:
            alt_n += 1
        elif not item["indels"]:
            ref_n += 1
        else:
            other += 1
    return ref_n, alt_n, other, True


def fractions(counts):
    ref_n, alt_n, other = counts
    exact = ref_n + alt_n
    total = exact + other
    return (
        alt_n / exact if exact else 0.0,
        other / total if total else 0.0,
        exact,
    )


result = run(["bcftools", "query", "-l", phased_vcf], capture_output=True)
samples = [line for line in result.stdout.splitlines() if line]
if len(samples) != 1:
    raise RuntimeError(f"Expected exactly one VCF sample, found {len(samples)}: {samples}")

full_records = []
for line in query_lines(
    full_vcf,
    "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER[\t%GT\t%DP\t%AD]\n",
):
    fields = line.split("\t")
    full_records.append({
        "chrom": fields[0],
        "pos": int(fields[1]),
        "ref": fields[2],
        "alt": fields[3],
        "qual": fields[4],
        "filter": fields[5],
        "original_gt": fields[6] if len(fields) > 6 else ".",
        "dp": fields[7] if len(fields) > 7 else ".",
        "ad": fields[8] if len(fields) > 8 else ".",
    })

phased = {}
for line in query_lines(
    phased_vcf,
    "%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%PS]\n",
):
    fields = line.split("\t")
    key = (fields[0], int(fields[1]), fields[2], fields[3])
    phased[key] = {
        "gt": fields[4] if len(fields) > 4 else ".",
        "ps": fields[5] if len(fields) > 5 else ".",
    }

biallelic_records = [r for r in full_records if "," not in r["alt"]]

with tempfile.TemporaryDirectory(prefix="hap_support_") as tmpdir:
    tmpdir = Path(tmpdir)
    sites_bed = tmpdir / "sites.bed"
    with sites_bed.open("w") as handle:
        for r in biallelic_records:
            handle.write(f'{r["chrom"]}\t{r["pos"] - 1}\t{r["pos"]}\n')

    pileups = {}
    for hp, bam in [("HP1", hp1_bam), ("HP2", hp2_bam)]:
        pileup_path = tmpdir / f"{hp}.pileup"
        with pileup_path.open("w") as out_handle:
            run([
                "samtools", "mpileup",
                "-B", "-Q", "0", "-q", "0",
                "-d", str(mpileup_max_depth),
                "-f", ref,
                "-l", str(sites_bed),
                bam,
            ], stdout=out_handle)
        data = {}
        with pileup_path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 5:
                    data[(fields[0], int(fields[1]))] = fields[4]
        pileups[hp] = data

    rows = []
    accepted_keys = set()
    unresolved_intervals = []

    for r in full_records:
        key = (r["chrom"], r["pos"], r["ref"], r["alt"])
        vtype = variant_type(r["ref"], r["alt"])
        phase = phased.get(key, {"gt": ".", "ps": "."})
        phased_gt = phase["gt"]
        ps = phase["ps"]

        hp_values = {}
        parser_supported = False
        if vtype != "MULTIALLELIC":
            for hp in ("HP1", "HP2"):
                bases = pileups[hp].get((r["chrom"], r["pos"]), "")
                ref_n, alt_n, other, supported = count_support(r["ref"], r["alt"], bases)
                parser_supported = parser_supported or supported
                alt_frac, other_frac, exact_depth = fractions((ref_n, alt_n, other))
                hp_values[hp] = {
                    "ref": ref_n,
                    "alt": alt_n,
                    "other": other,
                    "alt_frac": alt_frac,
                    "other_frac": other_frac,
                    "exact_depth": exact_depth,
                }
        else:
            for hp in ("HP1", "HP2"):
                hp_values[hp] = {
                    "ref": None, "alt": None, "other": None,
                    "alt_frac": None, "other_frac": None, "exact_depth": None,
                }

        status = "UNRESOLVED"
        reason = ""

        if vtype == "MULTIALLELIC":
            reason = "multiallelic_not_routinely_phased"
        elif key not in phased:
            reason = "not_present_in_phased_vcf"
        elif require_pass and r["filter"] != "PASS":
            reason = "non_pass_filter"
        elif not parser_supported or vtype in {"MNV", "COMPLEX"}:
            reason = "unsupported_complex_allele"
        else:
            h1 = hp_values["HP1"]
            h2 = hp_values["HP2"]
            enough_depth = (
                h1["exact_depth"] >= min_support_depth
                and h2["exact_depth"] >= min_support_depth
            )
            other_ok = (
                h1["other_frac"] <= max_other_fraction
                and h2["other_frac"] <= max_other_fraction
            )

            if phased_gt in {"1/1", "1|1"}:
                if (
                    enough_depth
                    and other_ok
                    and h1["alt_frac"] >= min_hom_alt_fraction
                    and h2["alt_frac"] >= min_hom_alt_fraction
                ):
                    status = "ACCEPT"
                    reason = "homozygous_alt_supported"
                else:
                    reason = "homozygous_alt_support_conflict"
            elif phased_gt in {"0|1", "1|0"}:
                expected = h2 if phased_gt == "0|1" else h1
                opposite = h1 if phased_gt == "0|1" else h2
                delta = expected["alt_frac"] - opposite["alt_frac"]
                if (
                    enough_depth
                    and other_ok
                    and expected["alt_frac"] >= min_het_alt_fraction
                    and delta >= min_het_delta
                ):
                    status = "ACCEPT"
                    reason = "phased_heterozygous_supported"
                else:
                    reason = "weak_or_conflicting_haplotype_support"
            elif phased_gt in {"0/1", "1/0"}:
                reason = "unphased_heterozygous"
            else:
                reason = f"unsupported_genotype_{phased_gt}"

        if status == "ACCEPT":
            accepted_keys.add(key)
        else:
            unresolved_intervals.append((r["chrom"], r["pos"] - 1, r["pos"] - 1 + len(r["ref"])))

        h1 = hp_values["HP1"]
        h2 = hp_values["HP2"]
        delta_abs = None
        if h1["alt_frac"] is not None and h2["alt_frac"] is not None:
            delta_abs = abs(h1["alt_frac"] - h2["alt_frac"])

        rows.append({
            **r,
            "type": vtype,
            "phased_gt": phased_gt,
            "ps": ps,
            "hp1": h1,
            "hp2": h2,
            "delta": delta_abs,
            "status": status,
            "reason": reason,
        })

    columns = [
        "CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "ORIGINAL_GT", "PHASED_GT", "PS", "TYPE",
        "HP1_REF", "HP1_ALT", "HP1_OTHER", "HP1_ALT_FRAC", "HP1_OTHER_FRAC",
        "HP2_REF", "HP2_ALT", "HP2_OTHER", "HP2_ALT_FRAC", "HP2_OTHER_FRAC",
        "ALT_FRAC_DELTA", "STATUS", "REASON",
    ]

    def fmt(value):
        if value is None:
            return "."
        if isinstance(value, float):
            return f"{value:.4f}"
        return str(value)

    def row_values(row):
        return [
            row["chrom"], row["pos"], row["ref"], row["alt"], row["qual"], row["filter"],
            row["original_gt"], row["phased_gt"], row["ps"], row["type"],
            row["hp1"]["ref"], row["hp1"]["alt"], row["hp1"]["other"],
            row["hp1"]["alt_frac"], row["hp1"]["other_frac"],
            row["hp2"]["ref"], row["hp2"]["alt"], row["hp2"]["other"],
            row["hp2"]["alt_frac"], row["hp2"]["other_frac"],
            row["delta"], row["status"], row["reason"],
        ]

    with support_tsv.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(fmt(v) for v in row_values(row)) + "\n")

    with uncertain_tsv.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            if row["status"] == "UNRESOLVED":
                handle.write("\t".join(fmt(v) for v in row_values(row)) + "\n")

    def merge_intervals(intervals):
        merged = []
        for chrom, start, end in sorted(intervals):
            if not merged or chrom != merged[-1][0] or start > merged[-1][2]:
                merged.append([chrom, start, end])
            else:
                merged[-1][2] = max(merged[-1][2], end)
        return merged

    with uncertain_bed.open("w") as handle:
        for chrom, start, end in merge_intervals(unresolved_intervals):
            handle.write(f"{chrom}\t{start}\t{end}\n")

    with accepted_bed.open("w") as handle:
        for chrom, pos, ref_allele, alt_allele in sorted(accepted_keys):
            handle.write(f"{chrom}\t{pos - 1}\t{pos}\t{ref_allele}\t{alt_allele}\n")

    plain_vcf = tmpdir / "consensus_ready.vcf"
    result = run(["bcftools", "view", "-Ov", phased_vcf], capture_output=True)
    with plain_vcf.open("w") as handle:
        for line in result.stdout.splitlines():
            if line.startswith("#"):
                handle.write(line + "\n")
                continue
            fields = line.split("\t")
            key = (fields[0], int(fields[1]), fields[3], fields[4])
            if key in accepted_keys:
                handle.write(line + "\n")

    run(["bcftools", "view", "-Oz", "-o", str(consensus_vcf), str(plain_vcf)])
    run(["tabix", "-f", "-p", "vcf", str(consensus_vcf)])

if not consensus_tbi.exists():
    raise RuntimeError(f"Expected tabix index was not created: {consensus_tbi}")
