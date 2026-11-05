rule minimap2:
    """
    Map reads to the reference genome with Minimap2
    """
    input:
        fastq= lambda wildcards: samples.loc[wildcards.sample, "fastq"],
        fasta="{wdir}/genome/{genome}.fasta",
        #html=lambda wildcards: run_qc_files(wildcards, "fastqc/{run}_sra_fastqc.html"),
        #qczip=lambda wildcards: run_qc_files(wildcards, "fastqc/{run}_sra_fastqc.zip"),
        #nanoplot=lambda wildcards: run_qc_files(
        #    wildcards, "nanoplot/{run}_NanoStats.txt"
        #),
        #longqc=lambda wildcards: run_qc_files(wildcards, "longqc/{run}"),
        #nanoplot_filtered="{wdir}/{sample}/nanoplot_filtered/{genome}_NanoStats.txt",
    output:
        sam=temp("{wdir}/{sample}/bam/{genome}_minimap2.sam"),
    conda:
        "../envs/minimap2.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_minimap2.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}.minimap2.tsv"
    params:
        minimap_ax=config["minimap_ax"],
        seed=config["seed"],
    shell:
        """
        minimap2 -ax {params.minimap_ax} --MD -2 \
        --seed {params.seed} --eqx \
        -t {resources.cpus_per_task} \
        -R "@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}" \
        --sam-hit-only {input.fasta} {input.fastq} > {output.sam}
        """



rule samtools_view:
    """
    Transform the sam file to a bam file

    One job per aligner, like rules samtools_sort and samtools_index below: a single job
    converting every alignment at once could not run when only one aligner is selected
    (config key `aligners`), because it would demand the other aligner's SAM.
    """
    input:
        sam="{wdir}/{sample}/bam/{genome}_{aligner}.sam",
    output:
        bam=temp("{wdir}/{sample}/bam/{genome}_{aligner}.bam"),
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}.samtools_view.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.samtools_view.tsv"
    shell:
        """
        samtools view -S -b {input.sam} > {output.bam}
        """


rule samtools_sort:
    """
    Sort the bam file
    """
    input:
        "{wdir}/{sample}/bam/{genome}_{aligner}.bam",
    output:
        "{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}.samtools_sort.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.samtools_sort.tsv"
    shell:
        """
        samtools sort {input} -o {output}
        """


rule samtools_index:
    """
    Create an index related file of the sorted bam file
    """
    input:
        "{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam",
    output:
        "{wdir}/{sample}/bam/{genome}_{aligner}_sorted.bam.bai",
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_{aligner}.samtools_index.log",
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_{aligner}.samtools_index.tsv"
    shell:
        """
        samtools index {input}
        """


#
# Alignment-level QC (samtools_stats, samtools_coverage) lives in
# rules/mapping_qc.smk, because it also applies to user-supplied BAM files.