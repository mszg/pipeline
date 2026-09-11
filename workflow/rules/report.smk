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
