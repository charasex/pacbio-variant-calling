import os
import re
import pandas as pd
from snakemake.exceptions import WorkflowError

configfile: "config/config.yaml"

# Import et validation des échantillons
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).fillna("")
if "fastq" not in samples.columns:
    raise WorkflowError("La colonne 'fastq' est absente de config/samples.tsv")

samples = samples.set_index("sample_name", drop=False).sort_index()

genomes = sorted(set(samples["genome"]))
if len(genomes) != 1:
    raise WorkflowError(f"Tous les échantillons doivent partager le même génome de référence. Trouvé : {genomes}")
genome = genomes[0]

individuals = sorted(set(samples["sample_name"]))
wdir = f"{config['datadir']}{genome}"

localrules: all

rule all:
    input:
        # 1. Contrôle qualité des reads bruts
        expand("{wdir}/{sample}/longqc", wdir=wdir, sample=individuals),
        # 2. Fichiers d'alignement BAM triés et indexés
        expand("{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai", wdir=wdir, sample=individuals, genome=genome),
        # 3. Contrôle qualité du mapping
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt", wdir=wdir, sample=individuals, genome=genome),
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv", wdir=wdir, sample=individuals, genome=genome)


rule longqc:
    # QC long read
    input:
        "{wdir}/read/{sample}.hifi_reads.bc2046.fastq.gz",
    output:
        directory("{wdir}/{sample}/longqc"),
    params:
        preset=config["longqc_preset"]
    threads: 8
    conda:
        "workflow/envs/longqc.yaml"
    log:
        "{wdir}/{sample}/logs/longqc.log"
    benchmark:
        "{wdir}/{sample}/benchmarks/longqc.tsv"
    shell:
        """
        python "$CONDA_PREFIX/share/LongQC/longQC.py" sampleqc \
        -x {params.preset} \
        -p $(( {threads} < 4 ? 4 : {threads} )) \
        -o {output} \
        --index 1G \
        {input} &> {log}
        """


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
        "workflow/envs/minimap2.yaml"
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
    # SAM vers BAM
    input:
        sam="{wdir}/{sample}/bam/{genome}_minimap2.sam"
    output:
        bam=temp("{wdir}/{sample}/bam/{genome}_minimap2.bam")
    conda:
        "workflow/envs/samtools.yaml"
    threads: 4
    shell:
        "samtools view -S -b -@ {threads} {input.sam} > {output.bam}"

rule samtools_sort:
    # Tri du fichier BAM par coordonnées génomiques
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    conda:
        "workflow/envs/samtools.yaml"
    threads: 8
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools_index:
    # Indexation du fichier BAM trié
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools index {input}"


rule samtools_stats:
    # Statistiques d'alignement globale
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools stats {input.bam} > {output}"

rule samtools_coverage:
    # Calcul de la profondeur de couverture par contig/chromosome
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools coverage {input.bam} -o {output}"