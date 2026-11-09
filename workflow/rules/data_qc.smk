rule nanoplot:
    input:
        fastq="{wdir}/read/{sample}.fastq.gz"
    output:
        "{wdir}/{sample}/nanoplot/{genome}_NanoStats.txt",
        "{wdir}/{sample}/nanoplot/{genome}_LengthvsQualityScatterPlot_dot.html",
        "{wdir}/{sample}/nanoplot/{genome}_LengthvsQualityScatterPlot_kde.html",
        "{wdir}/{sample}/nanoplot/{genome}_NanoPlot-report.html",
        "{wdir}/{sample}/nanoplot/{genome}_Non_weightedHistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{genome}_Non_weightedLogTransformed_HistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{genome}_WeightedHistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{genome}_WeightedLogTransformed_HistogramReadlength.html",
        "{wdir}/{sample}/nanoplot/{genome}_Yield_By_Length.html"
    threads: 8
    conda:
        "../envs/nanoplot.yaml"
    log:
        "{wdir}/{sample}/logs/{genome}_nanoplot.log"
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_nanoplot.tsv"
    shell:
        """
        NanoPlot --fastq {input.fastq} -t {threads} --tsv_stats \
        --outdir {wildcards.wdir}/{wildcards.sample}/nanoplot/ \
        --prefix '{wildcards.genome}_' --N50 --no_static --verbose \
        --title {wildcards.sample} &> {log}
        """