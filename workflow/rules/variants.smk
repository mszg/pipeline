rule clair3_call:
    input:
        bam="results/mapping/genes/{analysis}/{analysis}.merged.bam",
        bai="results/mapping/genes/{analysis}/{analysis}.merged.bam.bai",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai"
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
        model=lambda wc: config["clair3"]["model_path"],
        platform=lambda wc: config["clair3"]["platform"],
        extra=lambda wc: config["clair3"].get("extra", ""),
        outdir=lambda wc: f"results/variants/{wc.analysis}/clair3"
    shell:
        r"""
        mkdir -p {params.outdir:q} results/logs/variants
        {params.executable} \
          --bam_fn={input.bam:q} \
          --ref_fn={input.ref:q} \
          --threads={threads} \
          --platform={params.platform:q} \
          --model_path={params.model:q} \
          --output={params.outdir:q} \
          --sample_name={wildcards.analysis:q} \
          {params.extra} \
          > {log:q} 2>&1
        test -s {output.vcf:q}
        test -s {output.tbi:q}
        """


rule normalize_variants:
    input:
        vcf="results/variants/{analysis}/clair3/merge_output.vcf.gz",
        ref="results/reference/{analysis}/reference.fasta",
        fai="results/reference/{analysis}/reference.fasta.fai"
    output:
        vcf="results/variants/{analysis}/{analysis}.norm.vcf.gz",
        tbi="results/variants/{analysis}/{analysis}.norm.vcf.gz.tbi"
    conda:
        "../envs/variants.yaml"
    shell:
        r"""
        bcftools norm -f {input.ref:q} -m -any {input.vcf:q} -Oz -o {output.vcf:q}
        tabix -f -p vcf {output.vcf:q}
        """
