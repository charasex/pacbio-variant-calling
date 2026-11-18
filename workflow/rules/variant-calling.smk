rule bcftools_mpileup:
    """
    generate vcf with only variant site
    """
    input:
        ref="{wdir}/genome/{genome}.fasta",
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
            -X {params.configtools} -Oz -f {input.ref} {input.bam} | \
        bcftools call -mv -Oz -o {output.vcf}) 2> {log}
        """

rule bcftools_all_sites_roh:
    """
    Génère un VCF temporaire contenant les sites variants (0/1, 1/1) et invariants (0/0)
    bien couverts pour bcftools roh.
    """
    input:
        ref="{wdir}/genome/{primary}.fasta",
        bam="{wdir}/{sample}/bam/{primary}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{primary}_minimap2_sorted.bam.bai"
    output:
        vcf_all=temp("{wdir}/{sample}/mpileup/{primary}_bcftools_allsites.vcf.gz"),
        csi_all=temp("{wdir}/{sample}/mpileup/{primary}_bcftools_allsites.vcf.gz.csi")
    params:
        minMQ=config["minMQ"],
        minBQ=config["minBQ"],
        minQUAL=config["minQUAL"],
        configtools=config["configTools"],
        min_dp=config["DP"],
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{primary}_bcftools_allsites.log"
    threads: 8
    shell:
        """
        (bcftools mpileup --threads {threads} \
            -q {params.minMQ} -Q {params.minBQ} \
            -X {params.configtools} -Oz -f {input.ref} {input.bam} | \
        bcftools call -m -Oz | \
        bcftools filter -s LowQual -e 'QUAL<{params.minQUAL} || INFO/DP>{params.min_dp}' -Oz -o {output.vcf_all} && \
        bcftools index {output.vcf_all}) 2> {log}
        """



rule bfctools_index:
    input:
        vcf="{wdir}/{sample}/mpileup/{genome}_bcftools.vcf.gz"
    output:
        csi="{wdir}/{sample}/mpileup/{genome}_bcftools.vcf.gz.csi"
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_bcftools_index.log"
    threads: 1
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_bcftools_index.tsv"
    shell:
        """
        bcftools index {input.vcf} 2> {log}
        """
