"""Read-only report metrics from the pipeline's structured artifacts.

Missing files and malformed files are different evidence states. Counts are
observations, not molecules. Coverage is the union of intended reference bases.
"""
from collections import Counter, defaultdict
import csv
import gzip
import hashlib
import math
from pathlib import Path
import re
import statistics


class Evidence:
    def __init__(self):
        self.sources = []
        self.issues = []

    def read(self, path, parser, required=True, stage="analysis"):
        path = Path(path)
        source = {"file": str(path.resolve()), "stage": stage, "required": required}
        self.sources.append(source)
        if not path.exists():
            source["state"] = "MISSING"
            if required:
                self.issues.append({"status": "NOT ASSESSED", "reason": f"Missing {stage} evidence: {path}"})
            return None
        try:
            value = parser(path)
            source.update(state="READ", bytes=path.stat().st_size, mtime_ns=path.stat().st_mtime_ns,
                          sha256=hashlib.sha256(path.read_bytes()).hexdigest())
            return value
        except (ValueError, KeyError, TypeError, IndexError, OSError, EOFError, UnicodeError, csv.Error) as exc:
            source.update(state="PARSE ERROR", error=str(exc))
            self.issues.append({"status": "FAIL", "reason": f"Malformed {stage} evidence {path}: {exc}"})
            return None


def number(value, integer=False):
    if value in (None, "", ".", "NA", "N/A"):
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"Nonfinite numeric annotation {value!r}")
    if integer:
        if int(result) != result:
            raise ValueError(f"Expected integer, found {value!r}")
        return int(result)
    return result


def table(path, required=()):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not set(required) <= set(reader.fieldnames or []):
            raise ValueError(f"Expected columns {required}, found {reader.fieldnames}")
        if len(reader.fieldnames or []) != len(set(reader.fieldnames or [])):
            raise ValueError("Duplicate TSV column names")
        rows = list(reader)
    if any(None in row or any(v is None for v in row.values()) for row in rows):
        raise ValueError("Ragged TSV row")
    return rows


def metric_table(path):
    rows = table(path, ("metric", "value"))
    if len({r["metric"] for r in rows}) != len(rows):
        raise ValueError("Duplicate metric name")
    return {r["metric"]: r["value"] for r in rows}


def numeric_metrics(path, required):
    data = metric_table(path)
    if not set(required) <= data.keys():
        raise ValueError(f"Missing metrics {set(required) - data.keys()}")
    result = {k: number(data[k], True) for k in required}
    if any(v is None or v < 0 for v in result.values()):
        raise ValueError("Required counts must be nonnegative integers")
    return result


def primary_flagstat(path):
    text = path.read_text()
    result = {}
    for name in ("primary", "primary mapped", "secondary", "supplementary"):
        pattern = rf"^(\d+) \+ (\d+) {name}(?: \([^\n]*\))?$"
        found = re.findall(pattern, text, re.M)
        if len(found) != 1:
            raise ValueError(f"Missing or ambiguous flagstat {name}")
        result[name] = sum(map(int, found[0]))
    if result["secondary"] or result["supplementary"] or result["primary"] != result["primary mapped"]:
        raise ValueError("Cleaned BAM flagstat contains nonprimary or unmapped records")
    return result["primary mapped"]


def support_table(path):
    rows = table(path, ["CHROM", "POS", "REF", "ALT", "STATUS", "REASON"])
    keys = set()
    for row in rows:
        key = row["CHROM"], int(row["POS"]), row["REF"], row["ALT"]
        if key[1] < 1 or key in keys or row["STATUS"] not in {"ACCEPT", "UNRESOLVED"}:
            raise ValueError("Invalid, duplicate or unknown support record")
        keys.add(key)
    return rows


def hp_coverage_table(path):
    rows = table(path, ["metric", "HP1", "HP2"])
    if len({r["metric"] for r in rows}) != len(rows):
        raise ValueError("Duplicate haplotype coverage metric")
    required = {"target_positions", "callable_min_depth", "mean_depth", "min_depth", "callable_positions", "callable_fraction", "below_callable_depth"}
    if not required <= {r["metric"] for r in rows}:
        raise ValueError("Missing haplotype coverage metrics")
    counts = required - {"mean_depth", "callable_fraction"}
    parsed = {hp: {} for hp in ("HP1", "HP2")}
    for row in rows:
        for hp in ("HP1", "HP2"):
            value = number(row[hp], row["metric"] in counts)
            if value is None or value < 0:
                raise ValueError(f"Invalid {hp} coverage metric")
            parsed[hp][row["metric"]] = value
    for hp, values in parsed.items():
        total = values["target_positions"]
        if values["callable_min_depth"] <= 0 or values["min_depth"] > values["mean_depth"]:
            raise ValueError(f"Invalid {hp} depth threshold or minimum/mean ordering")
        if values["callable_positions"] + values["below_callable_depth"] != total:
            raise ValueError(f"{hp} callable/below-depth counts do not partition target positions")
        fraction = values["callable_fraction"]
        expected = values["callable_positions"] / total if total else 0
        if not 0 <= fraction <= 1 or abs(fraction - expected) > 0.000000500001:
            raise ValueError(f"Invalid {hp} callable fraction or disagreement with counts")
    return rows


def nano_reads(path):
    found = re.findall(r"^Number of reads:\s*([\d,.]+)\s*$", path.read_text(), re.M)
    if len(found) != 1:
        raise ValueError("Expected exactly one NanoStats Number of reads")
    value = number(found[0].replace(",", ""), True)
    if value is None or value < 0:
        raise ValueError("Invalid FASTQ count")
    return value


def merge_intervals(intervals):
    out = []
    for chrom, start, end in sorted(intervals):
        if start < 0 or end <= start:
            raise ValueError(f"Invalid interval {chrom}:{start}-{end}")
        if out and chrom == out[-1][0] and start <= out[-1][2]:
            out[-1] = (chrom, out[-1][1], max(end, out[-1][2]))
        else:
            out.append((chrom, start, end))
    return out


def manifest_targets(rows):
    intervals = []
    for row in rows:
        region = row.get("target_region", "")
        if not region:
            continue
        match = re.fullmatch(r"([^:\s]+):(\d+)-(\d+)", region)
        if not match:
            raise ValueError(f"Invalid manifest target {region!r}")
        chrom, start, end = match.groups()
        intervals.append((chrom, int(start) - 1, int(end)))
    return merge_intervals(intervals)


def bed(path):
    intervals = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split("\t")
        if len(fields) < 3:
            raise ValueError("BED needs at least three tab-separated fields")
        intervals.append((fields[0], int(fields[1]), int(fields[2])))
    return merge_intervals(intervals)


def positions(intervals):
    return {(chrom, p) for chrom, start, end in intervals for p in range(start, end)}


def regions(points):
    out = []
    for chrom, pos in sorted(points):
        if out and out[-1][0] == chrom and out[-1][2] == pos:
            out[-1][2] += 1
        else:
            out.append([chrom, pos, pos + 1])
    return out


def dense_depth(path, target):
    found = {}
    seen = set()
    for line in path.read_text().splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            raise ValueError("Expected samtools depth CHROM, POS, DEPTH")
        chrom, pos, depth = fields[0], int(fields[1]), int(fields[2])
        key = chrom, pos - 1
        if pos < 1 or depth < 0 or key in seen:
            raise ValueError("Invalid or duplicate depth position")
        seen.add(key)
        if key in target:
            found[key] = depth
    if found.keys() != target:
        raise ValueError(f"Dense -aa depth lacks {len(target - found.keys())} intended positions; missing is not zero")
    return found


def fasta_stats(path):
    lengths, ambiguities, n_counts, sequences = {}, {}, {}, {}
    name = None
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            name = line[1:].split()[0]
            if name in lengths:
                raise ValueError("Duplicate FASTA contig")
            lengths[name], ambiguities[name], n_counts[name] = 0, 0, 0
            sequences[name] = []
        elif line.strip():
            if name is None:
                raise ValueError("Sequence before FASTA header")
            seq = line.strip().upper()
            if not re.fullmatch(r"[ACGTRYSWKMBDHVN]+", seq):
                raise ValueError("Non-IUPAC FASTA bases")
            lengths[name] += len(seq)
            n_counts[name] += seq.count("N")
            ambiguities[name] += sum(b not in "ACGT" for b in seq)
            sequences[name].append(seq)
    if not lengths or any(length == 0 for length in lengths.values()):
        raise ValueError("Empty FASTA or contig without sequence")
    return {"lengths": lengths, "length": sum(lengths.values()), "n_bases": sum(n_counts.values()),
            "ambiguous_bases": sum(ambiguities.values()),
            "reference_ambiguities": {(c, i) for c, parts in sequences.items()
                                      for i, b in enumerate("".join(parts)) if b not in "ACGT"}}


def read_vcf(path, analysis):
    # Strict text parsing keeps FORMAT vectors and missing elements without
    # making DP=sum(AD) or imputing AF. Empty, header-valid VCFs are valid.
    opener = gzip.open if str(path).endswith(".gz") else open
    output, samples, keys = [], None, set()
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip().split("\t")[9:]
                if samples != [analysis]:
                    raise ValueError(f"VCF sample {samples!r} does not match {analysis!r}")
                continue
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if samples is None or len(f) != 10:
                raise ValueError("Missing VCF sample header or malformed record")
            pos = int(f[1])
            alts = f[4].split(",") if f[4] != "." else []
            key = f[0], pos, f[3], tuple(alts)
            if pos < 1 or key in keys:
                raise ValueError("Invalid/duplicate VCF record")
            keys.add(key)
            fmt = f[8].split(":")
            vals = f[9].split(":")
            if len(vals) != len(fmt) or len(set(fmt)) != len(fmt):
                raise ValueError("Malformed FORMAT/sample fields")
            data = dict(zip(fmt, vals))
            gt = data.get("GT", ".")
            gt_indices = [number(v, True) for v in re.split(r"[/|]", gt)]
            if any(v is not None and (v < 0 or v > len(alts)) for v in gt_indices):
                raise ValueError("GT allele outside REF/ALT catalogue")
            info = dict(field.split("=", 1) for field in f[7].split(";") if "=" in field)
            def vector(name):
                raw = data.get(name)
                return None if raw in (None, ".", "") else [number(x, name == "AD") for x in raw.split(",")]
            ad, af = vector("AD"), vector("AF")
            if any(v is not None and v < 0 for v in (ad or [])) or any(v is not None and not 0 <= v <= 1 for v in (af or [])):
                raise ValueError("Negative allele depths or AF outside [0,1]")
            if any(number(data.get(k)) is not None and number(data[k]) < 0 for k in ("DP",)):
                raise ValueError("Negative read depth")
            if ad is not None and len(ad) != len(alts) + 1:
                raise ValueError("AD must retain REF followed by all ALT alleles")
            if af is not None and len(af) != len(alts):
                raise ValueError("AF must retain an entry for each ALT allele")
            output.append({"contig": f[0], "position": pos, "REF": f[3], "ALT": f[4],
                           "genotype": None if gt == "." else gt, "FILTER": None if f[6] == "." else f[6],
                           "quality": number(f[5]), "depth": number(data.get("DP"), True),
                           "info_depth": number(info.get("DP"), True),
                           "allele_depths": ad, "allele_fractions": af,
                           "phase_set": None if data.get("PS") in (None, ".") else data["PS"],
                           "_key": key, "_alts": alts, "_gt": gt_indices,
                           "_heterozygous": len(gt_indices) == 2 and None not in gt_indices and gt_indices[0] != gt_indices[1],
                           "_phased": "|" in gt})
    if samples is None:
        raise ValueError("VCF has no #CHROM header")
    return output


def workflow_stages(config):
    w = config.get("workflow", {})
    v = bool(w.get("run_variant_calling", False))
    p = v and bool(w.get("run_phasing", False))
    h = p and bool(w.get("run_haplotype_reconstruction", False))
    return {"variants": v, "phasing": p, "haplotypes": h,
            "consensus": h and bool(w.get("run_consensus", False))}


def vcf_provenance(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    out = []
    with opener(path, "rt") as handle:
        for line in handle:
            if not line.startswith("##"):
                break
            key = line[2:].split("=", 1)[0].lower()
            if key in {"source", "reference", "cmdline", "commandline"} or key.endswith(("version", "command")):
                out.append(line.rstrip())
    return out


def collect_gene(rows, results, config, input_manifest=None):
    evidence = Evidence()
    gene, analysis = rows[0]["gene"], rows[0]["analysis"]
    units = [r["unit"] for r in rows]
    metrics = {"gene": gene, "analysis": analysis, "expected_amplicons": len(units),
               "observed_amplicons": 0, "_not_applicable": []}
    stages = workflow_stages(config)
    variant_rows, haplotype_rows = [], []
    if any(r.get("analysis_identity_ambiguous") for r in rows):
        evidence.issues.append({"status": "NOT ASSESSED", "reason": "Analysis output identity is shared by multiple manifest identities; pooled metrics are withheld"})
        return metrics, variant_rows, haplotype_rows, evidence
    if len(set(units)) != len(units):
        evidence.issues.append({"status": "FAIL", "reason": "Duplicate manifest sample/gene/amplicon unit; aggregation withheld"})
        return metrics, variant_rows, haplotype_rows, evidence
    reference = results / f"reference/{analysis}/reference.fasta"
    ref = evidence.read(reference, fasta_stats, stage="reference")
    try:
        intervals = manifest_targets(rows)
    except ValueError as exc:
        evidence.issues.append({"status": "FAIL", "reason": str(exc)})
        intervals = []
    target = positions(intervals)
    actual_bed = evidence.read(results / f"targets/{analysis}/{analysis}.bed", bed, stage="targets")
    if actual_bed is not None and actual_bed != intervals:
        evidence.issues.append({"status": "FAIL", "reason": "Output target BED disagrees with authoritative manifest union"})
        target = set()
    if ref and any(chrom not in ref["lengths"] or end > ref["lengths"].get(chrom, 0) for chrom, start, end in intervals):
        evidence.issues.append({"status": "FAIL", "reason": "Manifest target is outside staged reference"})
        target = set()
    if any(not row.get("target_region") for row in rows):
        evidence.issues.append({"status": "NOT ASSESSED", "reason": "One or more expected amplicons lack intended target coordinates; full gene target coverage is withheld"})
        target = set()
    if not target:
        evidence.issues.append({"status": "NOT ASSESSED", "reason": "Intended target intervals are unavailable or inconsistent; target coverage not assessed"})
    metrics.update(target_bases=len(target) if target else None,
                   reference_identifiers=", ".join(ref["lengths"]) if ref else None,
                   target_intervals="; ".join(f"{c}:{s+1}-{e}" for c,s,e in intervals), reference_file=str(reference.resolve()))
    input_counts, mapper_counts, mapped_counts, depths = [], [], [], []
    per_amplicon = []
    source_seen, overlap_sources = {}, []
    source_problem = False
    if input_manifest is not None:
        configured = {r["unit"]: r for r in rows}
        observed_units = set()
        for r in input_manifest:
            if r["unit"] not in units:
                continue
            observed_units.add(r["unit"])
            meta = configured[r["unit"]]
            source = (results.parent / r["source"]).resolve()
            requested = (results.parent / meta["fastq_input"]).resolve()
            is_fastq = str(requested).lower().endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))
            if any(r[key] != meta[key] for key in ("sample", "gene", "amplicon")) or (source != requested if is_fastq else not source.is_relative_to(requested)):
                source_problem = True
                evidence.issues.append({"status": "FAIL", "reason": f"Input provenance conflicts with current manifest for {r['unit']}: {source}"})
            if source in source_seen:
                overlap_sources.append(str(source))
            source_seen[source] = r["unit"]
        if observed_units != set(units):
            source_problem = True
            evidence.issues.append({"status": "NOT ASSESSED", "reason": "Input provenance lacks expected units: " + ", ".join(sorted(set(units)-observed_units))})
    elif len(units) > 1:
        evidence.issues.append({"status": "NOT ASSESSED", "reason": "Input manifest unavailable: cannot verify disjoint source-file aggregation"})
    for row in rows:
        unit = row["unit"]
        raw_count = evidence.read(results / f"qc/raw/{unit}/NanoStats.txt", nano_reads, stage=f"{unit} input QC")
        aqc = evidence.read(results / f"qc/alignment/{unit}/{unit}.alignment_qc.tsv",
                            lambda p: numeric_metrics(p, ["primary_records", "primary_mapped"]), stage=f"{unit} alignment QC")
        depth = evidence.read(results / f"mapping/amplicons/{unit}/{unit}.variant.depth.tsv",
                              lambda p: dense_depth(p, target), required=bool(target), stage=f"{unit} target depth") if target else None
        input_counts.append(raw_count)
        mapper_counts.append(aqc["primary_records"] if aqc else None)
        mapped_counts.append(aqc["primary_mapped"] if aqc else None)
        depths.append(depth)
        if aqc is not None:
            metrics["observed_amplicons"] += 1
            if aqc["primary_mapped"] > aqc["primary_records"]:
                evidence.issues.append({"status": "FAIL", "reason": f"{unit}: mapped count exceeds mapper-input primary count"})
        per_amplicon.append({"unit": unit, "amplicon": row["amplicon"], "barcode": row.get("barcode") or None,
                            "input_observations": raw_count, "alignment": aqc, "target_region": row.get("target_region"),
                            "phasing_min_length": row.get("phasing_min_length") or config.get("mapping", {}).get("phasing_min_length", 8000),
                            "core_depth_threshold": row.get("core_depth_threshold") or config.get("qc", {}).get("core_depth_threshold", 20)})
    metrics["_amplicons"] = per_amplicon
    metrics["source_file_count"] = len(source_seen) if input_manifest is not None else None
    if overlap_sources:
        evidence.issues.append({"status": "FAIL", "reason": "Repeated source FASTQ files across this gene's inputs; pooled read/depth metrics withheld: " + "; ".join(overlap_sources)})
    def total(values):
        return sum(values) if input_manifest is not None and not source_problem and not overlap_sources and all(v is not None for v in values) else None
    vc = evidence.read(results / f"mapping/genes/{analysis}/{analysis}.variant.flagstat.txt", primary_flagstat, stage="retained variant read counts")
    pc = evidence.read(results / f"mapping/genes/{analysis}/{analysis}.phasing.flagstat.txt", primary_flagstat, stage="retained phasing read counts")
    metrics.update(input_reads=total(input_counts), mapper_input_reads=total(mapper_counts), primary_mapped_reads=total(mapped_counts),
                   retained_variant_reads=vc if not overlap_sources and not source_problem and input_manifest is not None else None,
                   phasing_reads=pc if not overlap_sources and not source_problem and input_manifest is not None else None)
    numerator, denominator = metrics["primary_mapped_reads"], metrics["mapper_input_reads"]
    metrics["mapping_fraction"] = numerator / denominator if numerator is not None and denominator else None
    metrics["mapping_percent"] = metrics["mapping_fraction"] * 100 if metrics["mapping_fraction"] is not None else None
    metrics["mapping_denominator"] = "primary mapped / all primary alignment records (mapper-input observations; unmapped included; secondary/supplementary excluded)"
    metrics["_observations_definition"] = "FASTQ records across disjoint source files; BAM counts are primary observations, not deduplicated molecules; sample totals are not pooled across genes"
    rcfg = config.get("reporting", {})
    threshold = rcfg.get("assessment_depth", config.get("haplotypes", {}).get("callable_min_depth"))
    threshold = number(threshold, True)
    levels = set(rcfg.get("depth_thresholds", []))
    if threshold is not None:
        levels.add(threshold)
    if any(not isinstance(t, int) or isinstance(t, bool) or t <= 0 for t in levels):
        raise ValueError("reporting depth thresholds must be positive integers")
    metrics["assessment_depth"] = threshold
    if target and input_manifest is not None and not source_problem and not overlap_sources and all(d is not None for d in depths):
        summed = {pos: sum(d[pos] for d in depths) for pos in target}
        values = list(summed.values())
        metrics.update(mean_target_depth=sum(values)/len(values), median_target_depth=statistics.median(values), min_target_depth=min(values),
                       zero_coverage_bases=sum(v == 0 for v in values))
        metrics["_zero_coverage_regions_bed"] = regions(p for p, v in summed.items() if v == 0)
        metrics["zero_coverage_regions"] = len(metrics["_zero_coverage_regions_bed"])
        metrics["zero_coverage_intervals"] = "; ".join(f"{c}:{s+1}-{e}" for c,s,e in metrics["_zero_coverage_regions_bed"]) or "NONE (0 regions)"
        for t in sorted(levels):
            metrics[f"breadth_{t}x_fraction"] = sum(v >= t for v in values)/len(values)
        if threshold is not None:
            metrics["target_breadth_fraction"] = metrics[f"breadth_{threshold}x_fraction"]
            low = [p for p, v in summed.items() if v < threshold]
            metrics["low_coverage_bases"] = len(low)
            metrics["_low_coverage_regions_bed"] = regions(low)
            metrics["low_coverage_regions"] = len(metrics["_low_coverage_regions_bed"])
            metrics["low_coverage_intervals"] = "; ".join(f"{c}:{s+1}-{e}" for c,s,e in metrics["_low_coverage_regions_bed"]) or "NONE (0 regions)"
    norm_path = results / f"variants/{analysis}/{analysis}.norm.vcf.gz"
    norm = evidence.read(norm_path, lambda p: read_vcf(p, analysis), stage="variant calling") if stages["variants"] else None
    phase_path = results / f"phasing/{analysis}/{analysis}.phased.vcf.gz"
    phase = evidence.read(phase_path, lambda p: read_vcf(p, analysis), stage="phasing") if stages["phasing"] else None
    metrics.update(variant_records=len(norm) if norm is not None else None,
                   variants_file=str(norm_path.resolve()) if norm is not None else None,
                   phasing_file=str(phase_path.resolve()) if phase is not None else None)
    metrics["_vcf_provenance"] = []
    for artifact, parsed in ((norm_path, norm), (phase_path, phase)):
        if parsed is not None:
            metadata = evidence.read(artifact, vcf_provenance, stage="recorded VCF provenance")
            metrics["_vcf_provenance"].append({"file": str(artifact.resolve()), "headers": metadata})
    phase_by_key = {r["_key"]: r for r in phase} if phase is not None else {}
    if phase is not None and norm is not None:
        biallelic = {r["_key"]: r for r in norm if len(r["_alts"]) == 1}
        if biallelic.keys() != phase_by_key.keys():
            evidence.issues.append({"status": "FAIL", "reason": "Phased VCF catalogue differs from normalized biallelic records"})
        elif any(sorted(biallelic[k]["_gt"], key=str) != sorted(r["_gt"], key=str) for k,r in phase_by_key.items()):
            evidence.issues.append({"status": "FAIL", "reason": "Phasing changed an unphased genotype"})
    if norm is not None:
        metrics["multiallelic_records"] = sum(len(r["_alts"]) > 1 for r in norm)
        metrics["heterozygous_records"] = sum(r["_heterozygous"] for r in norm)
        metrics["nonbiallelic_heterozygous_records"] = sum(r["_heterozygous"] and len(r["_alts"]) != 1 for r in norm)
        metrics["pass_variant_records"] = sum(r["FILTER"] == "PASS" for r in norm)
        for r in norm:
            out = {k: v for k,v in r.items() if not k.startswith("_")}
            phased = phase_by_key.get(r["_key"])
            out.update(phased_genotype=phased["genotype"] if phased else None,
                       phase_set=phased["phase_set"] if phased else r["phase_set"],
                       allele_depths=",".join("." if v is None else str(v) for v in r["allele_depths"]) if r["allele_depths"] is not None else None,
                       allele_fractions=",".join("." if v is None else str(v) for v in r["allele_fractions"]) if r["allele_fractions"] is not None else None,
                       source_file=str(norm_path.resolve()))
            variant_rows.append(out)
    blocks = defaultdict(list)
    if phase is not None:
        het = [r for r in phase if r["_heterozygous"]]
        phased_hets = [r for r in het if r["_phased"] and r["phase_set"] is not None]
        if any(r["_phased"] and r["phase_set"] is None for r in het):
            evidence.issues.append({"status": "FAIL", "reason": "Phased heterozygote is missing its phase-set annotation"})
        metrics.update(eligible_heterozygous_records=len(het), phased_heterozygous_records=len(phased_hets),
                       unphased_heterozygous_records=len(het)-len(phased_hets), phased_fraction=len(phased_hets)/len(het) if het else None)
        for r in phased_hets:
            blocks[(r["contig"], r["phase_set"])].append(r["position"])
        metrics["phase_block_count"] = len(blocks)
        if norm is not None:
            metrics["all_unphased_heterozygous_records"] = metrics["heterozygous_records"] - len(phased_hets)
        metrics["phasing_applicability"] = "APPLICABLE" if het else "NOT APPLICABLE: no eligible heterozygous sites"
        if not het:
            metrics["_not_applicable"] += ["phased_fraction", "phase_block_count", "unphased_heterozygous_records"]
        for (chrom, ps), sites in sorted(blocks.items()):
            haplotype_rows.append({"row_type": "PHASE BLOCK", "contig": chrom, "phase_set": ps,
                                   "block_start": min(sites), "block_end": max(sites), "block_span_bp": max(sites)-min(sites)+1,
                                   "phased_heterozygous_records": len(sites), "source_file": str(phase_path.resolve()),
                                   "orientation": "HP labels are local to this gene/contig/phase set"})
    hapqc = results / f"qc/haplotypes/{analysis}"
    support = evidence.read(hapqc / f"{analysis}.variant_support.tsv", support_table, stage="allele support") if stages["haplotypes"] else None
    if support is not None:
        if norm is not None and {(r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]) for r in support} != {(r["contig"],r["position"],r["REF"],r["ALT"]) for r in norm}:
            evidence.issues.append({"status": "FAIL", "reason": "Allele-support catalogue differs from normalized VCF"})
        counts = Counter(r["STATUS"] for r in support)
        if set(counts) - {"ACCEPT", "UNRESOLVED"}:
            evidence.issues.append({"status": "FAIL", "reason": "Unknown support decision"})
        else:
            metrics.update(accepted_variants=counts["ACCEPT"], unresolved_variants=counts["UNRESOLVED"])
        lookup = {(r["CHROM"], int(r["POS"]), r["REF"], r["ALT"]): r for r in support}
        for row in variant_rows:
            s = lookup.get((row["contig"], row["position"], row["REF"], row["ALT"]), {})
            row.update(support_decision=s.get("STATUS"), support_reason=s.get("REASON"))
    tagging = evidence.read(hapqc / f"{analysis}.haplotag_qc.tsv", lambda p: numeric_metrics(p, ["total_reads", "hp1_reads", "hp2_reads", "assigned_reads", "unassigned_reads"]), stage="haplotagging") if stages["haplotypes"] else None
    if tagging is not None:
        if tagging["hp1_reads"] + tagging["hp2_reads"] != tagging["assigned_reads"] or tagging["assigned_reads"] + tagging["unassigned_reads"] != tagging["total_reads"]:
            evidence.issues.append({"status": "FAIL", "reason": "Haplotag counts do not partition input observations"})
        if metrics.get("phasing_reads") is not None and tagging["total_reads"] != metrics["phasing_reads"]:
            evidence.issues.append({"status": "FAIL", "reason": "Haplotag input denominator differs from phasing BAM count"})
        metrics.update(assigned_reads=tagging["assigned_reads"], unassigned_reads=tagging["unassigned_reads"],
                       assigned_read_fraction=tagging["assigned_reads"]/tagging["total_reads"] if tagging["total_reads"] else None,
                       haplotag_input_reads=tagging["total_reads"])
    hp_coverage = evidence.read(hapqc / f"{analysis}.coverage_qc.tsv", hp_coverage_table, stage="haplotype coverage") if stages["haplotypes"] else None
    hp_cov = {r["metric"]: r for r in hp_coverage} if hp_coverage else {}
    for hp in (1, 2):
        if not stages["haplotypes"] and not stages["consensus"]:
            continue
        h = {"row_type": "HAPLOTYPE", "haplotype": f"HP{hp}", "assigned_reads": tagging.get(f"hp{hp}_reads") if tagging else None,
             "phase_blocks": len(blocks) if phase is not None else None,
             "phased_heterozygous_records": metrics.get("phased_heterozygous_records"),
             "unphased_heterozygous_records": metrics.get("unphased_heterozygous_records"),
             "nonbiallelic_heterozygous_records": metrics.get("nonbiallelic_heterozygous_records"),
             "all_unphased_heterozygous_records": metrics.get("all_unphased_heterozygous_records"),
             "orientation": "Local gene/contig/phase-set labels; no linkage across genes or disconnected blocks inferred"}
        for name in ("target_positions", "callable_min_depth", "mean_depth", "min_depth", "callable_positions", "callable_fraction", "below_callable_depth"):
            try:
                h[name] = number(hp_cov[name][f"HP{hp}"]) if name in hp_cov else None
            except ValueError as exc:
                evidence.issues.append({"status": "FAIL", "reason": f"Malformed haplotype {name}: {exc}"})
                h[name] = None
        metrics[f"hp{hp}_min_depth"] = h.get("min_depth")
        if target and h.get("target_positions") is not None and h["target_positions"] != len(target):
            evidence.issues.append({"status": "FAIL", "reason": f"HP{hp} coverage target size disagrees with manifest union"})
        if stages["consensus"]:
            cons_path = results / f"consensus/{analysis}/{analysis}.haplotype{hp}.fasta"
            mask_path = results / f"consensus/{analysis}/{analysis}.HP{hp}.mask.bed"
            seq = evidence.read(cons_path, fasta_stats, stage=f"HP{hp} consensus")
            mask_intervals = evidence.read(mask_path, bed, stage=f"HP{hp} consensus mask")
            if ref and mask_intervals is not None and any(c not in ref["lengths"] or e > ref["lengths"].get(c, 0) for c,s,e in mask_intervals):
                evidence.issues.append({"status": "FAIL", "reason": f"HP{hp} mask is outside staged reference"})
                mask_intervals = None
            mask = positions(mask_intervals) if mask_intervals is not None else None
            if seq:
                h.update(consensus_length=seq["length"], fasta_n_bases=seq["n_bases"], fasta_ambiguous_bases=seq["ambiguous_bases"],
                         consensus_file=str(cons_path.resolve()), mask_file=str(mask_path.resolve()))
            if mask is not None and target:
                masked = len(mask & target)
                h.update(target_bases=len(target), target_masked_reference_bases=masked,
                         target_unmasked_reference_bases=len(target)-masked, target_unmasked_fraction=(len(target)-masked)/len(target),
                         outside_target_masked_bases=len(mask-target))
                metrics[f"hp{hp}_target_unmasked_fraction"] = h["target_unmasked_fraction"]
                metrics[f"hp{hp}_masked_target_bases"] = masked
                if ref:
                    h["unmasked_reference_ambiguous_target_bases"] = len((ref["reference_ambiguities"] & target) - mask)
            h["completeness_definition"] = "Fraction of intended reference-coordinate bases outside mask; supported deletions count as resolved. FASTA Ns include outside-target flanks. Not proof of biological phase."
        haplotype_rows.append(h)
    metrics["_stages"] = stages
    metrics["requested_stages"] = "; ".join(f"{name}: {'ENABLED' if on else 'DISABLED / NOT REQUESTED'}" for name,on in stages.items())
    metrics["_source_files"] = evidence.sources
    return metrics, variant_rows, haplotype_rows, evidence
