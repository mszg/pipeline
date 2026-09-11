rule amplicon_summary:
    input:
        nanostats=expand("results/qc/raw/{unit}/NanoStats.txt", unit=UNIT_IDS),
        flagstats=expand("results/mapping/amplicons/{unit}/{unit}.flagstat.txt", unit=UNIT_IDS),
        coverages=expand("results/mapping/amplicons/{unit}/{unit}.coverage.txt", unit=UNIT_IDS),
        manifest="results/summary/input_manifest.tsv"
    output:
        tsv="results/summary/amplicon_summary.tsv"
    params:
        samples_tsv=SAMPLES_TSV,
        units=",".join(UNIT_IDS)
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/make_amplicon_summary.py"


rule gene_summary:
    input:
        flagstats=expand("results/mapping/genes/{analysis}/{analysis}.flagstat.txt", analysis=ANALYSIS_IDS),
        coverages=expand("results/mapping/genes/{analysis}/{analysis}.coverage.txt", analysis=ANALYSIS_IDS),
        amplicon_summary="results/summary/amplicon_summary.tsv"
    output:
        tsv="results/summary/gene_summary.tsv"
    params:
        analyses=",".join(ANALYSIS_IDS)
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/make_gene_summary.py"
