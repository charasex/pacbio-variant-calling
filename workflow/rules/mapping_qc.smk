"""
Alignment-level QC.

These rules only need a sorted BAM and its index, so they run in both input modes:
after mapping (rules/mapping.smk) or on user-supplied BAMs (rules/bam_input.smk).
Their outputs feed rule final_report, which reads
mapping_QC/{genome}_{aligner}_mapping.stats.tsv for every selected aligner (config key
`aligners`).
"""


rule samtools_stats:
    """
    Mapping QC
    """
    input:
        bam="{wdir}/{genome}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{genome}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        stats="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_mapping.stats",
        stattsv="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_mapping.stats.tsv",
        plot="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_mapping_plot.html",
    conda:
        "../envs/samtools.yaml"
    log:
        "{wdir}/{genome}/logs/{aligner}.samtools_stats.log",
    benchmark:
        "{wdir}/{genome}/benchmarks/{aligner}.samtools_stats.tsv"
    params:
        wdir=wdir,
        genome=genome,
    shell:
        """
        mkdir -p {params.wdir}/{wildcards.genome}/mapping_QC
        samtools stats {input.bam} > {output.stats}
        cat {output.stats} | grep ^SN | cut -f 2- > {output.stattsv}
        # QC visualization
        plot-bamstats -p {params.wdir}/{wildcards.genome}/mapping_QC/{params.genome}_{wildcards.aligner}_mapping_plot {output.stats}
        """


rule samtools_coverage:
    """
    Mapping coverage along the genome
    """
    input:
        bam="{wdir}/{genome}/bam/{genome}_{aligner}_sorted.bam",
        bai="{wdir}/{genome}/bam/{genome}_{aligner}_sorted.bam.bai",
    output:
        coverage="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_coverage.tsv",
        coverage_hist="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_coverage_hist.txt",
        coverage_depth="{wdir}/{genome}/mapping_QC/{genome}_{aligner}_coverage_depthplot.txt",
    conda:
        "../envs/samtools_coverage.yaml"
    log:
        "{wdir}/{genome}/logs/{aligner}.samtools_coverage.log",
    benchmark:
        "{wdir}/{genome}/benchmarks/{aligner}.samtools_coverage.tsv"
    shell:
        """
        # store coverage in a tsv/tab-separated file
        samtools coverage {input.bam} > {output.coverage}

        # get coverage as a histogram
        samtools coverage {input.bam} --histogram > {output.coverage_hist}

        # get coverage as depth plot
        samtools coverage {input.bam} --plot-depth > {output.coverage_depth}
        """