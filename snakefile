import os
import re
import pandas as pd
from snakemake.exceptions import WorkflowError

configfile: "config/config.yaml"

# Import et validation des échantillons
samples = pd.read_csv(config["samples"], sep="\t", dtype=str).fillna("")
if "fastq" not in samples.columns:
    raise WorkflowError("the column 'fastq' is missing from config/samples.tsv")

samples = samples.set_index("sample_name", drop=False).sort_index()


genomes = sorted(set(samples["genome"]))
if len(genomes) != 1:
    raise WorkflowError(f"all samples must have the same reference genome.")
genome_name = genomes[0]

wdir = os.path.join(config["datadir"], genome_name)
individuals = sorted(set(samples["sample_name"]))

# def haplotype 
haplotypes = ["hap1", "hap2"]
primary_haplotype = config["primary_haplotype"]



localrules: all

rule all:
    input:
        # 1. control quality of reads
        expand("{wdir}/{sample}/nanoplot/{sample}_NanoPlot-report.html", wdir=wdir, sample=individuals),
        # 2. Bam alignements and index for hap1 AND hap2
        expand("{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai", wdir=wdir, sample=individuals, genome=haplotypes),
        # 3. control quality of mapping & coverage
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt", wdir=wdir, sample=individuals, genome=haplotypes),
        expand("{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv", wdir=wdir, sample=individuals, genome=haplotypes),
        expand("{wdir}/{sample}/mapping_QC/{genome}_mosdepth.mosdepth.summary.txt", wdir=wdir, sample=individuals, genome=haplotypes),
        expand("{wdir}/{sample}/mapping_QC/{genome}_mosdepth.regions.bed.gz", wdir=wdir, sample=individuals, genome=haplotypes),
        # 4. Variant calling uniquement sur l'haplotype principal
        expand("{wdir}/{sample}/mpileup/{primary}_bcftools.vcf.gz", wdir=wdir, sample=individuals, primary=primary_haplotype)

include: "workflow/rules/data_qc.smk"
include: "workflow/rules/mapping.smk"
include: "workflow/rules/mapping_qc.smk"
include: "workflow/rules/mappability.smk"
include: "workflow/rules/variant-calling.smk"