rule amplicon_summary:
    input:
        nanostats=expand("results/qc/raw/{unit}/NanoStats.txt", unit=UNIT_IDS),
        flagstats=expand("results/mapping/amplicons/{unit}/{unit}.flagstat.txt", unit=UNIT_IDS),
        raw_coverages=expand("results/mapping/amplicons/{unit}/{unit}.coverage.txt", unit=UNIT_IDS),
        alignment_qc=expand("results/qc/alignment/{unit}/{unit}.alignment_qc.tsv", unit=UNIT_IDS),
        analysis_coverages=expand("results/mapping/amplicons/{unit}/{unit}.variant.coverage.txt", unit=UNIT_IDS),
        core_intervals=expand("results/mapping/amplicons/{unit}/{unit}.core_intervals.tsv", unit=UNIT_IDS),
        manifest="results/summary/input_manifest.tsv"
    output:
        tsv="results/summary/amplicon_summary.tsv"
    params:
        samples_tsv=SAMPLES_TSV,
        units=",".join(UNIT_IDS),
        default_phasing_min_length=int(config.get("mapping", {}).get("phasing_min_length", 8000)),
        default_core_depth_threshold=int(config.get("qc", {}).get("core_depth_threshold", 20))
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/make_amplicon_summary.py"


rule gene_summary:
    input:
        raw_flagstats=expand("results/mapping/genes/{analysis}/{analysis}.flagstat.txt", analysis=ANALYSIS_IDS),
        raw_coverages=expand("results/mapping/genes/{analysis}/{analysis}.coverage.txt", analysis=ANALYSIS_IDS),
        variant_flagstats=expand("results/mapping/genes/{analysis}/{analysis}.variant.flagstat.txt", analysis=ANALYSIS_IDS),
        variant_coverages=expand("results/mapping/genes/{analysis}/{analysis}.variant.coverage.txt", analysis=ANALYSIS_IDS),
        phasing_flagstats=expand("results/mapping/genes/{analysis}/{analysis}.phasing.flagstat.txt", analysis=ANALYSIS_IDS),
        amplicon_summary="results/summary/amplicon_summary.tsv"
    output:
        tsv="results/summary/gene_summary.tsv"
    params:
        analyses=",".join(ANALYSIS_IDS)
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/make_gene_summary.py"


def evaluation_report_inputs(wc):
    files = list(analysis_final_targets) + [SAMPLES_TSV]
    files += [str(Path(workflow.basedir) / "workflow/scripts" / name) for name in
              ("evaluation_metrics.py", "evaluation_policy.py", "evaluation_workbook.py")]
    for analysis in ANALYSIS_IDS:
        files += [f"results/reference/{analysis}/reference.fasta",
                  f"results/mapping/genes/{analysis}/{analysis}.variant.flagstat.txt",
                  f"results/mapping/genes/{analysis}/{analysis}.phasing.flagstat.txt"]
        if run_variants:
            files.append(f"results/variants/{analysis}/{analysis}.norm.vcf.gz")
        if run_phasing:
            files.append(f"results/phasing/{analysis}/{analysis}.phased.vcf.gz")
        if run_consensus:
            files += [f"results/consensus/{analysis}/{analysis}.HP{hp}.mask.bed" for hp in (1, 2)]
    files += [f"results/mapping/amplicons/{unit}/{unit}.variant.depth.tsv" for unit in UNIT_IDS]
    return list(dict.fromkeys(files))


rule evaluation_report:
    input:
        evaluation_report_inputs
    output:
        xlsx="results/reports/pipeline_evaluation.xlsx",
        json="results/reports/pipeline_evaluation.json"
    params:
        effective_config=lambda wc: dict(config),
        repository=str(Path(workflow.basedir)),
        script_dir=str(Path(workflow.basedir) / "workflow/scripts")
    conda:
        "../envs/evaluation_report.yaml"
    script:
        "../scripts/evaluation_report.py"
