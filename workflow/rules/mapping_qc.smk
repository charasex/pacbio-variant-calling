"""
Alignment-level QC.

These rules only need a sorted BAM and its index, so they run in both input modes:
after mapping (rules/mapping.smk).
Their outputs feed rule final_report, which reads
"""

rule samtools_stats:
    # stat
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_stats.txt"
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_samtools_stats.log"
    shell:
        "samtools stats {input.bam} > {output} 2> {log}"

rule samtools_coverage:
    # coverage of each contig/chromosome
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        "{wdir}/{sample}/mapping_QC/{genome}_minimap2_coverage.tsv"
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_samtools_coverage.log"
    shell:
        "samtools coverage {input.bam} -o {output} 2> {log}"

