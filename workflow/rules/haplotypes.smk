rule haplotag_reads:
    input:
        vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi",
        bam="results/mapping/genes/{analysis}/{analysis}.phasing.bam",
        bai="results/mapping/genes/{analysis}/{analysis}.phasing.bam.bai",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai"
    output:
        bam="results/haplotypes/{analysis}/{analysis}.haplotagged.bam",
        bai="results/haplotypes/{analysis}/{analysis}.haplotagged.bam.bai"
    log:
        "results/logs/haplotypes/{analysis}.haplotag.log"
    conda:
        "../envs/haplotypes.yaml"
    threads:
        lambda wc: int(config["haplotypes"]["threads"])
    shell:
        r"""
        mkdir -p results/haplotypes/{wildcards.analysis} results/logs/haplotypes
        whatshap haplotag \
          --reference {input.ref:q} \
          --output {output.bam:q} \
          {input.vcf:q} {input.bam:q} \
          > {log:q} 2>&1
        samtools index -@ {threads} {output.bam:q}
        """


rule split_haplotagged_reads:
    input:
        bam="results/haplotypes/{analysis}/{analysis}.haplotagged.bam",
        bai="results/haplotypes/{analysis}/{analysis}.haplotagged.bam.bai"
    output:
        hp1_bam="results/haplotypes/{analysis}/{analysis}.HP1.bam",
        hp1_bai="results/haplotypes/{analysis}/{analysis}.HP1.bam.bai",
        hp2_bam="results/haplotypes/{analysis}/{analysis}.HP2.bam",
        hp2_bai="results/haplotypes/{analysis}/{analysis}.HP2.bam.bai"
    conda:
        "../envs/haplotypes.yaml"
    threads:
        lambda wc: int(config["haplotypes"]["threads"])
    shell:
        r"""
        samtools view -@ {threads} -b -d HP:1 {input.bam:q} -o {output.hp1_bam:q}
        samtools view -@ {threads} -b -d HP:2 {input.bam:q} -o {output.hp2_bam:q}
        samtools index -@ {threads} {output.hp1_bam:q}
        samtools index -@ {threads} {output.hp2_bam:q}
        """


rule haplotag_qc:
    input:
        bam="results/haplotypes/{analysis}/{analysis}.haplotagged.bam",
        bai="results/haplotypes/{analysis}/{analysis}.haplotagged.bam.bai"
    output:
        tsv="results/qc/haplotypes/{analysis}/{analysis}.haplotag_qc.tsv"
    conda:
        "../envs/haplotypes.yaml"
    script:
        "../scripts/haplotag_qc.py"


rule haplotype_coverage_qc:
    input:
        hp1_bam="results/haplotypes/{analysis}/{analysis}.HP1.bam",
        hp1_bai="results/haplotypes/{analysis}/{analysis}.HP1.bam.bai",
        hp2_bam="results/haplotypes/{analysis}/{analysis}.HP2.bam",
        hp2_bai="results/haplotypes/{analysis}/{analysis}.HP2.bam.bai",
        bed="results/targets/{analysis}/{analysis}.bed"
    output:
        summary="results/qc/haplotypes/{analysis}/{analysis}.coverage_qc.tsv",
        hp1_low="results/qc/haplotypes/{analysis}/{analysis}.HP1.low_depth.bed",
        hp2_low="results/qc/haplotypes/{analysis}/{analysis}.HP2.low_depth.bed"
    params:
        min_depth=lambda wc: int(config["haplotypes"]["callable_min_depth"])
    conda:
        "../envs/haplotypes.yaml"
    script:
        "../scripts/haplotype_coverage_qc.py"


rule haplotype_variant_support:
    input:
        full_vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        full_tbi="results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi",
        phased_vcf="results/phasing/{analysis}/{analysis}.phased.vcf.gz",
        phased_tbi="results/phasing/{analysis}/{analysis}.phased.vcf.gz.tbi",
        hp1_bam="results/haplotypes/{analysis}/{analysis}.HP1.bam",
        hp1_bai="results/haplotypes/{analysis}/{analysis}.HP1.bam.bai",
        hp2_bam="results/haplotypes/{analysis}/{analysis}.HP2.bam",
        hp2_bai="results/haplotypes/{analysis}/{analysis}.HP2.bam.bai",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai"
    output:
        support="results/qc/haplotypes/{analysis}/{analysis}.variant_support.tsv",
        uncertain="results/qc/haplotypes/{analysis}/{analysis}.uncertain_variants.tsv",
        uncertain_bed="results/qc/haplotypes/{analysis}/{analysis}.uncertain_regions.bed",
        accepted_bed="results/qc/haplotypes/{analysis}/{analysis}.accepted_variants.bed",
        consensus_vcf="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz",
        consensus_tbi="results/variants/{analysis}/{analysis}.consensus_ready.vcf.gz.tbi"
    params:
        min_support_depth=lambda wc: int(config["haplotypes"]["min_support_depth"]),
        min_het_alt_fraction=lambda wc: float(config["haplotypes"]["min_het_alt_fraction"]),
        min_het_delta=lambda wc: float(config["haplotypes"]["min_het_delta"]),
        min_hom_alt_fraction=lambda wc: float(config["haplotypes"]["min_hom_alt_fraction"]),
        max_other_fraction=lambda wc: float(config["haplotypes"]["max_other_fraction"]),
        require_pass=lambda wc: bool(config["haplotypes"].get("require_pass", True)),
        mpileup_max_depth=lambda wc: int(config["haplotypes"].get("mpileup_max_depth", 100000))
    conda:
        "../envs/haplotypes.yaml"
    script:
        "../scripts/haplotype_variant_support.py"
