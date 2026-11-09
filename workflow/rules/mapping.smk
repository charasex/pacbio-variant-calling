rule minimap2:
    # Alignement des reads 
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample, "fastq"],
        fasta=config["ref"],
    output:
        sam=temp("{wdir}/{sample}/bam/{genome}_minimap2.sam")
    params:
        minimap_ax=config["minimap_ax"],
        seed=config["seed"]
    threads: 16
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_minimap2.log"
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.minimap2.tsv"
    shell:
        """
        minimap2 -ax {params.minimap_ax} --MD -2 \
        --seed {params.seed} --eqx \
        -t {threads} \
        -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
        --sam-hit-only {input.fasta} {input.fastq} > {output.sam} 2> {log}
        """

rule samtools_view:
    # SAM --> BAM
    input:
        sam="{wdir}/{sample}/bam/{genome}_minimap2.sam"
    output:
        bam=temp("{wdir}/{sample}/bam/{genome}_minimap2.bam")
    conda:
        "../envs/samtools.yaml"
    threads: 4
    log:
        "{wdir}/{sample}/logs/{genome}_samtools_view.log"
    shell:
        "samtools view -S -b -@ {threads} {input.sam} > {output.bam} 2> {log}"

rule samtools_sort:
    # sort of BAM file 
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    conda:
        "../envs/samtools.yaml"
    threads: 10
    log:
        "{wdir}/{sample}/logs/{genome}_samtools_sort.log"
    shell:
        "samtools sort -@ {threads} {input} -o {output} 2> {log}"

rule samtools_index:
    # index of bam file sorted
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_samtools_index.log"
    shell:
        "samtools index {input} 2> {log}"

#
# Alignment-level QC (samtools_stats, samtools_coverage) lives in
# rules/mapping_qc.smk, because it also applies to user-supplied BAM files.