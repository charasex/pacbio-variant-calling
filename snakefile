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
        # 1. quality control of reads
        expand("{wdir}/{sample}/nanoplot/{genome}_NanoPlot-report.html", wdir=wdir, sample=individuals, genome=genome),
        # 2. Bam file of alignment sorted and indexed
        expand("{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai", wdir=wdir, sample=individuals, genome=genome),
        # 3. quality control of mapping
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt", wdir=wdir, sample=individuals, genome=genome),
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv", wdir=wdir, sample=individuals, genome=genome),
        # for variant calling
        expand("{wdir}/{sample}/mpileup/{genome}_bcftools.vcf.gz", wdir=wdir, sample=individuals, genome=genome)


include: "workflow/rules/data_qc.smk"
include: "workflow/rules/mapping.smk"
include: "workflow/rules/mapping_qc.smk"
include: "workflow/rules/variant-calling.smk"