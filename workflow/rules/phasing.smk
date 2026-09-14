rule phase_variants:
    input:
        vcf="results/variants/{analysis}/{analysis}.phasing_ready.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.phasing_ready.vcf.gz.tbi",
        bam="results/mapping/genes/{analysis}/{analysis}.phasing.bam",
        bai="results/mapping/genes/{analysis}/{analysis}.phasing.bam.bai",
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


rule phasing_qc:
    input:
        full_vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        phasing_ready_vcf="results/variants/{analysis}/{analysis}.phasing_ready.vcf.gz",
        phased_vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        phased_tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi"
    output:
        tsv="results/qc/phasing/{analysis}/{analysis}.phasing_qc.tsv"
    conda:
        "../envs/phasing.yaml"
    script:
        "../scripts/phasing_qc.py"
