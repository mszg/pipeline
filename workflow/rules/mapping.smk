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


rule amplicon_alignment_qc:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        summary="results/qc/alignment/{unit}/{unit}.alignment_qc.tsv",
        length_bins="results/qc/alignment/{unit}/{unit}.length_mapq.tsv"
    conda:
        "../envs/bam_qc.yaml"
    script:
        "../scripts/alignment_qc.py"


# Small-variant analysis BAM: primary mapped reads only, MAPQ-filtered.
# Target coordinates are NOT applied at BAM level in v0.4; they are passed
# explicitly to Clair3 as a BED calling interval so the BAM retains alignment context.
rule filter_variant_bam:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        bam="results/mapping/amplicons/{unit}/{unit}.variant.bam"
    log:
        "results/logs/mapping/{unit}.variant_filter.log"
    params:
        min_mapq=lambda wc: int(config["mapping"].get("variant_min_mapq", 30)),
        min_read_length=0,
        target_region="",
        allow_empty=False
    conda:
        "../envs/bam_qc.yaml"
    script:
        "../scripts/filter_bam.py"


rule index_variant_bam:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.variant.bam"
    output:
        bai="results/mapping/amplicons/{unit}/{unit}.variant.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


rule variant_bam_coverage:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.variant.bam"
    output:
        txt="results/mapping/amplicons/{unit}/{unit}.variant.coverage.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools coverage {input.bam:q} > {output.txt:q}"


rule variant_bam_depth:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.variant.bam"
    output:
        tsv="results/mapping/amplicons/{unit}/{unit}.variant.depth.tsv"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools depth -aa {input.bam:q} > {output.tsv:q}"


rule core_intervals:
    input:
        depth="results/mapping/amplicons/{unit}/{unit}.variant.depth.tsv"
    output:
        tsv="results/mapping/amplicons/{unit}/{unit}.core_intervals.tsv"
    params:
        threshold=core_depth_threshold_for_unit
    conda:
        "../envs/report.yaml"
    script:
        "../scripts/core_intervals.py"


# Phasing BAM: same primary/MAPQ cleanup plus a configurable long-read threshold.
# It also retains alignment context outside the calling BED; WhatsHap receives only
# variants produced inside the explicit Clair3 target intervals.
rule filter_phasing_bam:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.sorted.bam"
    output:
        bam="results/mapping/amplicons/{unit}/{unit}.phasing.bam"
    log:
        "results/logs/mapping/{unit}.phasing_filter.log"
    params:
        min_mapq=lambda wc: int(config["mapping"].get("phasing_min_mapq", 30)),
        min_read_length=phasing_min_length_for_unit,
        target_region="",
        allow_empty=True
    conda:
        "../envs/bam_qc.yaml"
    script:
        "../scripts/filter_bam.py"


rule index_phasing_bam:
    input:
        bam="results/mapping/amplicons/{unit}/{unit}.phasing.bam"
    output:
        bai="results/mapping/amplicons/{unit}/{unit}.phasing.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


# Keep a raw merged gene BAM for traceability/QC.
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


# This is the BAM used for Clair3.
rule merge_gene_variant_bams:
    input:
        bams=variant_bams_for_analysis
    output:
        bam="results/mapping/genes/{analysis}/{analysis}.variant.bam"
    log:
        "results/logs/mapping/{analysis}.variant_merge.log"
    conda:
        "../envs/mapping.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p results/mapping/genes/{wildcards.analysis} results/logs/mapping
        samtools merge -@ {threads} -f {output.bam:q} {input.bams:q} 2> {log:q}
        """


rule index_gene_variant_bam:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.variant.bam"
    output:
        bai="results/mapping/genes/{analysis}/{analysis}.variant.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


rule gene_variant_flagstat:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.variant.bam"
    output:
        txt="results/mapping/genes/{analysis}/{analysis}.variant.flagstat.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam:q} > {output.txt:q}"


rule gene_variant_coverage:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.variant.bam"
    output:
        txt="results/mapping/genes/{analysis}/{analysis}.variant.coverage.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools coverage {input.bam:q} > {output.txt:q}"


# This long-read-enriched BAM is used for WhatsHap.
rule merge_gene_phasing_bams:
    input:
        bams=phasing_bams_for_analysis
    output:
        bam="results/mapping/genes/{analysis}/{analysis}.phasing.bam"
    log:
        "results/logs/mapping/{analysis}.phasing_merge.log"
    conda:
        "../envs/mapping.yaml"
    threads: 4
    shell:
        r"""
        mkdir -p results/mapping/genes/{wildcards.analysis} results/logs/mapping
        samtools merge -@ {threads} -f {output.bam:q} {input.bams:q} 2> {log:q}
        """


rule index_gene_phasing_bam:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.phasing.bam"
    output:
        bai="results/mapping/genes/{analysis}/{analysis}.phasing.bam.bai"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools index {input.bam:q}"


rule gene_phasing_flagstat:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.phasing.bam"
    output:
        txt="results/mapping/genes/{analysis}/{analysis}.phasing.flagstat.txt"
    conda:
        "../envs/mapping.yaml"
    shell:
        "samtools flagstat {input.bam:q} > {output.txt:q}"
