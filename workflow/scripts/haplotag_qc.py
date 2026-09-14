import subprocess
from collections import Counter
from pathlib import Path

bam = str(snakemake.input.bam)
out = Path(str(snakemake.output.tsv))
out.parent.mkdir(parents=True, exist_ok=True)

hp_counts = Counter()
ps_counts = Counter()
total = 0
assigned = 0

proc = subprocess.Popen(
    ["samtools", "view", bam],
    stdout=subprocess.PIPE,
    text=True,
)
assert proc.stdout is not None
for line in proc.stdout:
    fields = line.rstrip("\n").split("\t")
    total += 1
    hp = None
    ps = None
    for tag in fields[11:]:
        if tag.startswith("HP:i:"):
            hp = tag.split(":", 2)[2]
        elif tag.startswith("PS:i:"):
            ps = tag.split(":", 2)[2]
    if hp in {"1", "2"}:
        hp_counts[hp] += 1
        assigned += 1
        if ps is not None:
            ps_counts[ps] += 1

return_code = proc.wait()
if return_code != 0:
    raise subprocess.CalledProcessError(return_code, ["samtools", "view", bam])

unassigned = total - assigned
assigned_fraction = assigned / total if total else 0.0
hp1_fraction = hp_counts["1"] / assigned if assigned else 0.0
hp2_fraction = hp_counts["2"] / assigned if assigned else 0.0

rows = [
    ("total_reads", total),
    ("hp1_reads", hp_counts["1"]),
    ("hp2_reads", hp_counts["2"]),
    ("assigned_reads", assigned),
    ("unassigned_reads", unassigned),
    ("assigned_fraction", f"{assigned_fraction:.6f}"),
    ("hp1_fraction_of_assigned", f"{hp1_fraction:.6f}"),
    ("hp2_fraction_of_assigned", f"{hp2_fraction:.6f}"),
    ("number_of_phase_sets", len(ps_counts)),
]

if ps_counts:
    dominant_ps, dominant_n = ps_counts.most_common(1)[0]
    rows.extend([
        ("dominant_phase_set", dominant_ps),
        ("dominant_phase_set_reads", dominant_n),
    ])
    for ps, count in sorted(ps_counts.items(), key=lambda x: (int(x[0]) if x[0].isdigit() else x[0])):
        rows.append((f"phase_set_{ps}_reads", count))
else:
    rows.extend([
        ("dominant_phase_set", "."),
        ("dominant_phase_set_reads", 0),
    ])

with out.open("w") as handle:
    handle.write("metric\tvalue\n")
    for metric, value in rows:
        handle.write(f"{metric}\t{value}\n")
