rule clair3_call:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.variant.bam",
        bai="results/mapping/genes/{analysis}/{analysis}.variant.bam.bai",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai",
        bed="results/targets/{analysis}/{analysis}.bed"
    output:
        vcf="results/variants/{analysis}/clair3/merge_output.vcf.gz",
        tbi="results/variants/{analysis}/clair3/merge_output.vcf.gz.tbi"
    log:
        "results/logs/variants/{analysis}.clair3.log"
    conda:
        "../envs/clair3.yaml"
    threads:
        lambda wc: int(config["clair3"]["threads"])
    params:
        executable=lambda wc: config["clair3"]["executable"],
        model_name=lambda wc: config["clair3"]["model_name"],
        platform=lambda wc: config["clair3"]["platform"],
        extra=lambda wc: config["clair3"].get("extra", ""),
        bed_arg=clair3_bed_arg,
        outdir=lambda wc: f"results/variants/{wc.analysis}/clair3"
    shell:
        r"""
        mkdir -p {params.outdir:q} results/logs/variants
        {params.executable} \
          --bam_fn={input.bam:q} \
          --ref_fn={input.ref:q} \
          --threads={threads} \
          --platform={params.platform:q} \
          --model_path="$(dirname "$(command -v run_clair3.sh)")/models/{params.model_name}" \
          --output={params.outdir:q} \
          --sample_name={wildcards.analysis:q} \
          {params.bed_arg} \
          {params.extra} \
          > {log:q} 2>&1
        test -s {output.vcf:q}
        test -s {output.tbi:q}
        """


rule normalize_variants:
    input:
        vcf="results/variants/{analysis}/clair3/merge_output.vcf.gz",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai",
        bed="results/targets/{analysis}/{analysis}.bed"
    output:
        vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi"
    params:
        target_arg=bcftools_target_arg
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        # Preserve multiallelic records in the complete normalized catalogue.
        # Splitting a 1/2 site into duplicate biallelic positions can make one
        # record invisible to WhatsHap, as observed during ABO validation.
        bcftools norm -f {input.ref:q} {input.vcf:q} -Ou \
          | bcftools view {params.target_arg} -Oz -o {output.vcf:q}
        tabix -f -p vcf {output.vcf:q}
        """


rule make_phasing_ready_variants:
    input:
        vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi"
    output:
        vcf="results/variants/{analysis}/{analysis}.phasing_ready.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.phasing_ready.vcf.gz.tbi"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        # WhatsHap receives only records with exactly one ALT allele.
        # The complete normalized VCF remains available for reporting so that
        # multiallelic candidates are preserved rather than silently lost.
        bcftools view -m2 -M2 -Oz -o {output.vcf:q} {input.vcf:q}
        tabix -f -p vcf {output.vcf:q}
        """
