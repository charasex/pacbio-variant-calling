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
        # 1. QC of reads 
        expand("{wdir}/{genome}/longqc/QC_result.html", wdir=wdir, sample=individuals, genome=genome),
        # 2. BAM
        expand("{wdir}/{genome}/bam/{sample}_minimap2_sorted.bam.bai", wdir=wdir, sample=individuals, genome=genome),
        # 3. QC of mapping
        expand("{wdir}/{genome}/mapping_QC/{sample}_minimap2_stats.txt", wdir=wdir, sample=individuals, genome=genome),
        expand("{wdir}/{genome}/mapping_QC/{sample}_minimap2_coverage.tsv", wdir=wdir, sample=individuals, genome=genome)