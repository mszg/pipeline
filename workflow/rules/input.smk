rule combine_fastq_chunks:
    input:
        fastqs=unit_fastqs
    output:
        fastq="results/input/{unit}/{unit}.combined.fastq.gz"
    log:
        "results/logs/input/{unit}.combine.log"
    script:
        "../scripts/combine_fastqs.py"


rule input_manifest:
    input:
        fastqs=lambda wc: [
            path
            for unit in UNIT_IDS
            for path in discover_fastqs(UNITS[unit]["fastq_input"])
        ]
    output:
        tsv="results/summary/input_manifest.tsv"
    run:
        import csv
        from pathlib import Path

        Path(output.tsv).parent.mkdir(parents=True, exist_ok=True)
        with open(output.tsv, "w", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t")
            writer.writerow(["sample", "gene", "amplicon", "unit", "chunk", "source"])
            for unit in UNIT_IDS:
                row = UNITS[unit]
                for i, path in enumerate(discover_fastqs(row["fastq_input"]), start=1):
                    writer.writerow([row["sample"], row["gene"], row["amplicon"], unit, i, path])
