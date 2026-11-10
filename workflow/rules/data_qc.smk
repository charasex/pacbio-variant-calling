rule nanoplot:
    input:
        fastq=lambda wildcards: samples.loc[wildcards.sample, "fastq"]
    output:
        "{wdir}/{sample}/nanoplot/{sample}_NanoStats.txt",
        "{wdir}/{sample}/nanoplot/{sample}_LengthvsQualityScatterPlot_dot.html",
        "{wdir}/{sample}/nanoplot/{sample}_LengthvsQualityScatterPlot_kde.html",
        "{wdir}/{sample}/nanoplot/{sample}_NanoPlot-report.html",
        "{wdir}/{sample}/nanoplot/{sample}_Non_weightedHistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{sample}_Non_weightedLogTransformed_HistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{sample}_WeightedHistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{sample}_WeightedLogTransformed_HistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{sample}_Yield_By_Length.html"
    threads: 8
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/{sample}/logs/{sample}_nanoplot.log"
    benchmark:
        "{wdir}/{sample}/benchmarks/{sample}_nanoplot.tsv"
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {threads} --tsv_stats \
        --outdir {wildcards.wdir}/{wildcards.sample}/nanoplot/ \
        --prefix '{wildcards.sample}_' --N50 --no_static --verbose \
        --title {wildcards.sample} &> {log}
        """