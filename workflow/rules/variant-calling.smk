rule bcftools_mpileup:
    input:
        ref=config["ref"],
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        vcf="{wdir}/{sample}/mpileup/{genome}_bcftools.vcf.gz"
    params:
        minMQ=config["minMQ"],
        minBQ=config["minBQ"],
        configtools=config["configTools"]
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_bcftools.log"
    threads: 8
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_bcftools.tsv"
    shell:
        """
        (bcftools mpileup --threads {threads} \
        -q {params.minMQ} -Q {params.minBQ} \
        -X {params.configtools} -Oz -f {input.ref} {input.bam} |\
        bcftools call -mv -Oz -o {output.vcf}) 2> {log}
        """
