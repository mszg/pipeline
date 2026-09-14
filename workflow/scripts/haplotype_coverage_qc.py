import subprocess
from collections import defaultdict
from pathlib import Path

bed_path = Path(str(snakemake.input.bed))
hp_bams = {
    "HP1": str(snakemake.input.hp1_bam),
    "HP2": str(snakemake.input.hp2_bam),
}
threshold = int(snakemake.params.min_depth)
summary_path = Path(str(snakemake.output.summary))
low_paths = {
    "HP1": Path(str(snakemake.output.hp1_low)),
    "HP2": Path(str(snakemake.output.hp2_low)),
}


def read_bed(path):
    intervals = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            chrom, start, end, *_ = line.rstrip("\n").split("\t")
            intervals.append((chrom, int(start), int(end)))
    return intervals


def merge_intervals(intervals):
    by_chrom = defaultdict(list)
    for chrom, start, end in intervals:
        by_chrom[chrom].append((start, end))
    merged = []
    for chrom, values in by_chrom.items():
        current = None
        for start, end in sorted(values):
            if current is None or start > current[1]:
                if current is not None:
                    merged.append((chrom, current[0], current[1]))
                current = [start, end]
            else:
                current[1] = max(current[1], end)
        if current is not None:
            merged.append((chrom, current[0], current[1]))
    return merged


def write_low_depth_bed(path, positions):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        current = None
        for chrom, pos1 in positions:
            start = pos1 - 1
            end = pos1
            if current is None:
                current = [chrom, start, end]
            elif chrom == current[0] and start == current[2]:
                current[2] = end
            else:
                handle.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")
                current = [chrom, start, end]
        if current is not None:
            handle.write(f"{current[0]}\t{current[1]}\t{current[2]}\n")


targets = merge_intervals(read_bed(bed_path))
target_positions = sum(end - start for _, start, end in targets)
results = {}

for hp, bam in hp_bams.items():
    depth_map = {}
    proc = subprocess.Popen(
        ["samtools", "depth", "-aa", "-b", str(bed_path), bam],
        stdout=subprocess.PIPE,
        text=True,
    )
    assert proc.stdout is not None
    for line in proc.stdout:
        chrom, pos, depth = line.rstrip("\n").split("\t")[:3]
        depth_map[(chrom, int(pos))] = int(depth)
    return_code = proc.wait()
    if return_code != 0:
        raise subprocess.CalledProcessError(return_code, ["samtools", "depth", bam])

    values = []
    low_positions = []
    for chrom, start, end in targets:
        for pos1 in range(start + 1, end + 1):
            depth = depth_map.get((chrom, pos1), 0)
            values.append(depth)
            if depth < threshold:
                low_positions.append((chrom, pos1))

    if len(values) != target_positions:
        raise RuntimeError(
            f"Internal coverage accounting error for {hp}: {len(values)} != {target_positions}"
        )

    write_low_depth_bed(low_paths[hp], low_positions)
    results[hp] = {
        "mean_depth": sum(values) / len(values) if values else 0.0,
        "min_depth": min(values) if values else 0,
        "max_depth": max(values) if values else 0,
        "zero_depth": sum(d == 0 for d in values),
        "depth_lt_10": sum(d < 10 for d in values),
        "depth_lt_20": sum(d < 20 for d in values),
        "depth_lt_50": sum(d < 50 for d in values),
        "depth_lt_100": sum(d < 100 for d in values),
        "below_callable_depth": len(low_positions),
        "callable_positions": target_positions - len(low_positions),
        "callable_fraction": (target_positions - len(low_positions)) / target_positions if target_positions else 0.0,
    }

summary_path.parent.mkdir(parents=True, exist_ok=True)
with summary_path.open("w") as handle:
    handle.write("metric\tHP1\tHP2\n")
    handle.write(f"target_positions\t{target_positions}\t{target_positions}\n")
    handle.write(f"callable_min_depth\t{threshold}\t{threshold}\n")
    for metric in [
        "mean_depth", "min_depth", "max_depth", "zero_depth",
        "depth_lt_10", "depth_lt_20", "depth_lt_50", "depth_lt_100",
        "below_callable_depth", "callable_positions", "callable_fraction",
    ]:
        a = results["HP1"][metric]
        b = results["HP2"][metric]
        if isinstance(a, float):
            a = f"{a:.6f}"
        if isinstance(b, float):
            b = f"{b:.6f}"
        handle.write(f"{metric}\t{a}\t{b}\n")
