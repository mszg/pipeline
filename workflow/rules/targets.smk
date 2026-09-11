rule build_analysis_target_bed:
    input:
        samples=SAMPLES_TSV
    output:
        bed="results/targets/{analysis}/{analysis}.bed",
        manifest="results/targets/{analysis}/{analysis}.target_regions.tsv"
    params:
        analysis=lambda wc: wc.analysis
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/make_target_bed.py"
