rule stage_reference:
    input:
        ref=analysis_reference
    output:
        ref="results/reference/{analysis}/reference.fasta"
    shell:
        r"""
        mkdir -p results/reference/{wildcards.analysis}
        cp {input.ref:q} {output.ref:q}
        """


rule faidx_reference:
    input:
        ref="results/reference/{analysis}/reference.fasta"
    output:
        fai="results/reference/{analysis}/reference.fasta.fai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools faidx {input.ref:q}"


rule minimap2_index:
    input:
        ref="results/reference/{analysis}/reference.fasta"
    output:
        mmi="results/reference/{analysis}/reference.mmi"
    conda:
        "../envs/mapping.yaml"
    shell:
        "minimap2 -d {output.mmi:q} {input.ref:q}"


rule align_amplicon:
    input:
        reads=reads_for_mapping,
        ref=staged_reference_for_unit,
        mmi=staged_mmi_for_unit
    output:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    log:
        "results/logs/mapping/{unit}.minimap2.log"
    conda:
        "../envs/mapping.yaml"
    threads:
        lambda wc: int(config["mapping"]["threads"])
    params:
        preset=lambda wc: config["mapping"]["preset"],
        rg=lambda wc: (
            f"@RG\\tID:{wc.unit}\\tSM:{analysis_for_unit(wc.unit)}\\tPL:ONT"
        )
    shell:
        r"""
        mkdir -p results/mapping/amplicons/{wildcards.unit} results/logs/mapping
        minimap2 -t {threads} -ax {params.preset} -R {params.rg:q} \
          {input.mmi:q} {input.reads:q} 2> {log:q} \
          | samtools sort -@ {threads} -o {output.bam:q} -
        """


rule index_amplicon_bam:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        bai="results/mapping/amplicons/{unit}/{unit}.sorted.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


rule amplicon_flagstat:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        txt="results/mapping/amplicons/{unit}/{unit}.flagstat.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam:q} > {output.txt:q}"


rule amplicon_coverage:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        txt="results/mapping/amplicons/{unit}/{unit}.coverage.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools coverage {input.bam:q} > {output.txt:q}"


rule merge_gene_bams:
    input:
        bams=amplicon_bams_for_analysis
    output:
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam"
    log:
        "results/logs/mapping/{analysis}.merge.log"
    conda:
        "../envs/mapping.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p results/mapping/genes/{wildcards.analysis} results/logs/mapping
        samtools merge -@ {threads} -f {output.bam:q} {input.bams:q} 2> {log:q}
        """


rule index_gene_bam:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam"
    output:
        bai="results/mapping/genes/{analysis}/{analysis}.merged.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


rule gene_flagstat:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam"
    output:
        txt="results/mapping/genes/{analysis}/{analysis}.flagstat.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam:q} > {output.txt:q}"


rule gene_coverage:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam"
    output:
        txt="results/mapping/genes/{analysis}/{analysis}.coverage.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools coverage {input.bam:q} > {output.txt:q}"
