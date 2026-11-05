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
        expand("{wdir}/{sample}/longqc/QC_result.html", wdir=wdir, sample=individuals),
        # 2. Fichiers d'alignement BAM triés et indexés
        expand("{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai", wdir=wdir, sample=individuals, genome=genome),
        # 3. Contrôle qualité du mapping
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt", wdir=wdir, sample=individuals, genome=genome),
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv", wdir=wdir, sample=individuals, genome=genome)

# ------------------------------------------------------------------------------
# ÉTAPE 1 : QC DES READS BRUTS (LongQC)
# ------------------------------------------------------------------------------
rule longqc:
    """
    Contrôle qualité des reads longs bruts avec LongQC
    """
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample, "fastq"]
    output:
        html="{wdir}/{sample}/longqc/QC_result.html",
        directory("{wdir}/{sample}/longqc"),
    params:
        preset=config["longqc_preset"]
    threads: 8
    conda:
        "envs/longqc.yaml"
    log:
        "{wdir}/{sample}/logs/longqc.log"
    shell:
        """
        longQC.py sample -x {params.preset} -o {output.outdir} -p {threads} {input.fastq} > {log} 2>&1
        """

# ------------------------------------------------------------------------------
# ÉTAPE 2 : ALIGNEMENT / MAPPING (Minimap2 + Samtools)
# ------------------------------------------------------------------------------
rule minimap2:
    """
    Alignement des reads locaux avec Minimap2 (-R épuré de read_group_platform)
    """
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample, "fastq"],
        fasta=config["ref"]
    output:
        sam=temp("{wdir}/{sample}/bam/{genome}_minimap2.sam")
    params:
        minimap_ax=config["minimap_ax"],
        seed=config["seed"]
    threads: 16
    conda:
        "envs/minimap2.yaml"
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
    """
    Conversion SAM vers BAM temporaire
    """
    input:
        sam="{wdir}/{sample}/bam/{genome}_minimap2.sam"
    output:
        bam=temp("{wdir}/{sample}/bam/{genome}_minimap2.bam")
    conda:
        "envs/samtools.yaml"
    threads: 4
    shell:
        "samtools view -S -b -@ {threads} {input.sam} > {output.bam}"

rule samtools_sort:
    """
    Tri du fichier BAM par coordonnées génomiques
    """
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    conda:
        "envs/samtools.yaml"
    threads: 8
    shell:
        "samtools sort -@ {threads} {input} -o {output}"

rule samtools_index:
    """
    Indexation du fichier BAM trié
    """
    input:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam"
    output:
        "{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools index {input}"

# ------------------------------------------------------------------------------
# ÉTAPE 3 : QC DU MAPPING
# ------------------------------------------------------------------------------
rule samtools_stats:
    """
    Statistiques d'alignement globale
    """
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools stats {input.bam} > {output}"

rule samtools_coverage:
    """
    Calcul de la profondeur de couverture par contig/chromosome
    """
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv"
    conda:
        "envs/samtools.yaml"
    shell:
        "samtools coverage {input.bam} -o {output}"