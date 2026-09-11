from pathlib import Path
import re


def parse_nanostats(path):
    out = {"mean_read_length": "NA", "mean_read_quality": "NA", "n50": "NA", "reads": "NA"}
    text = Path(path).read_text(errors="replace")
    patterns = {
        "mean_read_length": r"Mean read length:\s*([0-9.,]+)",
        "mean_read_quality": r"Mean read quality:\s*([0-9.,]+)",
        "n50": r"Read length N50:\s*([0-9.,]+)",
        "reads": r"Number of reads:\s*([0-9.,]+)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text, flags=re.I)
        if m:
            out[key] = m.group(1).replace(",", "")
    return out


def parse_flagstat(path):
    text = Path(path).read_text(errors="replace")
    total = "NA"
    mapped_count = "NA"
    mapped_percent = "NA"
    m = re.search(r"^(\d+) \+ \d+ in total", text, flags=re.M)
    if m:
        total = m.group(1)
    m = re.search(r"^(\d+) \+ \d+ mapped \(([0-9.]+)%", text, flags=re.M)
    if m:
        mapped_count = m.group(1)
        mapped_percent = m.group(2)
    return total, mapped_count, mapped_percent


def parse_coverage(path):
    rows = []
    for line in Path(path).read_text(errors="replace").splitlines():
        if not line or line.startswith("#"):
            continue
        vals = line.split("\t")
        if len(vals) < 7:
            continue
        start = int(vals[1])
        end = int(vals[2])
        length = max(0, end - start + 1)
        rows.append((length, float(vals[5]), float(vals[6])))

    if not rows:
        return "NA", "NA"

    total_len = sum(x[0] for x in rows)
    if total_len == 0:
        return "NA", "NA"
    coverage = sum(length * cov for length, cov, depth in rows) / total_len
    depth = sum(length * depth for length, cov, depth in rows) / total_len
    return f"{coverage:.4f}", f"{depth:.4f}"
