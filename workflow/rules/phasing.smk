rule phase_variants:
    input:
        vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi",
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam",
        bai="results/mapping/genes/{analysis}/{analysis}.merged.bam.bai",
        ref="results/reference/{analysis}/reference.fasta"
    output:
        vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi"
    log:
        "results/logs/phasing/{analysis}.whatshap.log"
    conda:
        "../envs/phasing.yaml"
    threads:
        lambda wc: int(config["phasing"]["threads"])
    shell:
        r"""
        mkdir -p results/phasing/{wildcards.analysis} results/logs/phasing
        whatshap phase \
          --reference {input.ref:q} \
          --output {output.vcf:q} \
          {input.vcf:q} {input.bam:q} \
          > {log:q} 2>&1
        tabix -f -p vcf {output.vcf:q}
        """
