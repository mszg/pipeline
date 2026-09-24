rule consensus_mask_hp1:
    input:
        fai="results/reference/{analysis}/reference.fasta.fai",
        target_bed="results/targets/{analysis}/{analysis}.bed",
        uncertain_bed="results/qc/haplotypes/{analysis}/{analysis}.uncertain_regions.bed",
        support_tsv="results/qc/haplotypes/{analysis}/{analysis}.variant_support.tsv",
        low_depth_bed="results/qc/haplotypes/{analysis}/{analysis}.HP1.low_depth.bed"
    output:
        bed="results/consensus/{analysis}/{analysis}.HP1.mask.bed"
    conda:
        "../envs/haplotypes.yaml"
    params:
        haplotype=1,
        callable_min_depth=lambda wc: int(config["haplotypes"]["callable_min_depth"])
    script:
        "../scripts/make_consensus_mask.py"


rule consensus_mask_hp2:
    input:
        fai="results/reference/{analysis}/reference.fasta.fai",
        target_bed="results/targets/{analysis}/{analysis}.bed",
        uncertain_bed="results/qc/haplotypes/{analysis}/{analysis}.uncertain_regions.bed",
        support_tsv="results/qc/haplotypes/{analysis}/{analysis}.variant_support.tsv",
        low_depth_bed="results/qc/haplotypes/{analysis}/{analysis}.HP2.low_depth.bed"
    output:
        bed="results/consensus/{analysis}/{analysis}.HP2.mask.bed"
    conda:
        "../envs/haplotypes.yaml"
    params:
        haplotype=2,
        callable_min_depth=lambda wc: int(config["haplotypes"]["callable_min_depth"])
    script:
        "../scripts/make_consensus_mask.py"


rule consensus_haplotype1:
    input:
        vcf="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz.tbi",
        ref="results/reference/{analysis}/reference.fasta",
        mask="results/consensus/{analysis}/{analysis}.HP1.mask.bed"
    output:
        fasta="results/consensus/{analysis}/{analysis}.haplotype1.fasta"
    conda:
        "../envs/haplotypes.yaml"
    shell:
        r"""
        mkdir -p results/consensus/{wildcards.analysis}
        if [ -s {input.mask:q} ]; then
            bcftools consensus -f {input.ref:q} -H 1 -m {input.mask:q} {input.vcf:q} > {output.fasta:q}
        else
            bcftools consensus -f {input.ref:q} -H 1 {input.vcf:q} > {output.fasta:q}
        fi
        """


rule consensus_haplotype2:
    input:
        vcf="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz.tbi",
        ref="results/reference/{analysis}/reference.fasta",
        mask="results/consensus/{analysis}/{analysis}.HP2.mask.bed"
    output:
        fasta="results/consensus/{analysis}/{analysis}.haplotype2.fasta"
    conda:
        "../envs/haplotypes.yaml"
    shell:
        r"""
        mkdir -p results/consensus/{wildcards.analysis}
        if [ -s {input.mask:q} ]; then
            bcftools consensus -f {input.ref:q} -H 2 -m {input.mask:q} {input.vcf:q} > {output.fasta:q}
        else
            bcftools consensus -f {input.ref:q} -H 2 {input.vcf:q} > {output.fasta:q}
        fi
        """
