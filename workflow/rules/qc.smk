rule nanoplot_raw:
    input:
        fastq="results/input/{unit}/{unit}.combined.fastq.gz"
    output:
        stats="results/qc/raw/{unit}/NanoStats.txt"
    log:
        "results/logs/qc/{unit}.nanoplot.log"
    conda:
        "../envs/qc.yaml"
    threads: 4
    params:
        static_option=lambda wc: "--no_static" if config.get("qc", {}).get("nanoplot_no_static", False) else ""
    shell:
        r"""
        mkdir -p results/qc/raw/{wildcards.unit} results/logs/qc
        NanoPlot --fastq {input.fastq:q} \
            --outdir results/qc/raw/{wildcards.unit} \
            --threads {threads} \
            {params.static_option} \
            > {log:q} 2>&1
        test -s {output.stats:q}
        """


rule filter_reads:
    input:
        fastq="results/input/{unit}/{unit}.combined.fastq.gz"
    output:
        fastq="results/filtered/{unit}/{unit}.filtered.fastq.gz"
    log:
        "results/logs/qc/{unit}.filtlong.log"
    conda:
        "../envs/qc.yaml"
    params:
        min_length=lambda wc: config["filtering"]["min_length"],
        max_length=lambda wc: config["filtering"]["max_length"],
        min_mean_q=lambda wc: config["filtering"]["min_mean_q"],
        keep_percent=lambda wc: config["filtering"]["keep_percent"]
    shell:
        r"""
        mkdir -p results/filtered/{wildcards.unit} results/logs/qc
        filtlong \
            --min_length {params.min_length} \
            --max_length {params.max_length} \
            --min_mean_q {params.min_mean_q} \
            --keep_percent {params.keep_percent} \
            {input.fastq:q} 2> {log:q} | gzip -c > {output.fastq:q}
        """
