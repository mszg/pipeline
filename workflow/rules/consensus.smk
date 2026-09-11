rule consensus_haplotype1:
    input:
        vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi",
        ref="results/reference/{analysis}/reference.fasta"
    output:
        fasta="results/consensus/{analysis}/{analysis}.haplotype1.fasta"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        mkdir -p results/consensus/{wildcards.analysis}
        bcftools consensus -f {input.ref:q} -H 1 {input.vcf:q} > {output.fasta:q}
        """


rule consensus_haplotype2:
    input:
        vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi",
        ref="results/reference/{analysis}/reference.fasta"
    output:
        fasta="results/consensus/{analysis}/{analysis}.haplotype2.fasta"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        mkdir -p results/consensus/{wildcards.analysis}
        bcftools consensus -f {input.ref:q} -H 2 {input.vcf:q} > {output.fasta:q}
        """
