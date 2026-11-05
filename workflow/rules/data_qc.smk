rule longqc:
    """
    Long-read quality control with LongQC, one SRA ?????? run at a time.

    LongQC: Fukasawa et al. (2020) G3: Genes, Genomes, Genetics 10(4):1193-1196.
    sampleqc gives read length, quality, GC content, low-complexity masking, coverage
    from read overlaps and adapter analysis at the read ends.
    -x      Platform preset (config key longqc_preset; from sequencing_technology:
            hifi -> pb-hifi, ont -> ont-ligation)
    -p      Number of CPUs. LongQC exits below 4, so it is raised to 4
    -o      Output folder. LongQC exits if it exists, so it is removed first

    LongQC analyses a sample of 5000 reads (-n default), drawn with a fixed seed (7),
    so the result is reproducible. Adapter trimming (-c) is not used: rules
    hifiadapterfilt and porechop_abi remove adapters.

    The output is the whole folder: web_summary.html loads its figures from figs/.

    longQC.py is run from $CONDA_PREFIX/share/LongQC, not from the copy bioconda
    puts in bin/. LongQC resolves its sdust and minimap2-coverage binaries relative
    to its own source file, a layout the bioconda packages do not provide, so
    envs/longqc.post-deploy.sh builds them into that tree. Running the copy in bin/
    makes every sdust call fail silently and the run dies ~30 min later on an empty
    longqc_sdust.txt. The script has the details. It installs longQC.py without a
    shebang, so it runs with python.

    Knwon issue:
    https://github.com/yfukasawa/LongQC/issues/28
    https://github.com/yfukasawa/LongQC/issues/76
    the error means that the process ran out of memory. This can be addressed by lowering the index parameter (-i/--index) for the minimap step. For my dataset lowering it from the default '4G' to '1G' (python3 longQC.py sampleqc -x pb-sequel -i 1G ... ) solved the issue.
    """
    input:
        "{wdir}/read/{sample}.hifi_reads.bc2046.fastq.gz",
    output:
        directory("{wdir}/{sample}/longqc"),
    threads: workflow.cores
    conda:
        "../envs/longqc.yaml"
    log:
        "{wdir}/{genome}/logs/longqc.log",
    benchmark:
        "{wdir}/{genome}/benchmarks/longqc.tsv"
    params:
        longqc_preset=config["longqc_preset"],
    shell:
        """
        python "$CONDA_PREFIX/share/LongQC/longQC.py" sampleqc \
        -x {params.longqc_preset} \
        -p $(( {threads} < 4 ? 4 : {threads} )) \
        -o {output} \
        --index 1G \
        {input} &> {log}
        """


rule filter_reads_chopper:
    """
    Filter long reads with chopper
    --headcrop      Trim N nucleotides from the start of a read
    --maxlength     Sets a maximum read length
    -l, --minlength     Sets a minimum read length
    -q, --quality       Sets a minimum Phred average quality score
    --tailcrop      Trim N nucleotides from the end of a read
    """
    input:
        reads="{wdir}/{genome}/read/{sample}.hifi_reads.bc2046.fastq.gz",
    output:
        filtered_reads=temp("{wdir}/{genome}/read/{sample}_filtered.fastq.gz"),
    conda:
        "../envs/chopper.yaml"
    log:
        "{wdir}/{genome}/logs/filter_reads_chopper.log",
    benchmark:
        "{wdir}/{genome}/benchmarks/filter_reads_chopper.tsv"
    params:
        chopper_quality=config["chopper_quality"],
        chopper_minlength=config["chopper_minlength"],
        chopper_maxlength=config["chopper_maxlength"],
        chopper_headcrop=config["chopper_headcrop"],
        chopper_tailcrop=config["chopper_tailcrop"],
    shell:
        """
        chopper -q {params.chopper_quality} \
        -l {params.chopper_minlength} \
        --maxlength {params.chopper_maxlength} \
        --headcrop {params.chopper_headcrop} \
        --tailcrop {params.chopper_tailcrop} \
        --threads {resources.cpus_per_task} \
        -i {input.reads} | gzip > {output.filtered_reads}
        """



# rule nanoplot_after_filtering:
#     """
#     Quality control after filtering long reads
#     """
#     input:
#         fastq="{wdir}/{sample}/fastq/{genome}_filtered.fastq.gz",
#     output:
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_NanoStats.txt",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_dot.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_LengthvsQualityScatterPlot_kde.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_NanoPlot-report.html",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedHistogramReadlength.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_Non_weightedLogTransformed_HistogramReadlength.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedHistogramReadlength.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_WeightedLogTransformed_HistogramReadlength.png",
#         "{wdir}/{sample}/nanoplot_filtered/{genome}_Yield_By_Length.html",
#         # "{wdir}/{sample}/nanoplot_filtered/{genome}_Yield_By_Length.png"
#     conda:
#         "../envs/nanoplot.yaml"
#     log:
#         "{wdir}/{sample}/logs/{genome}_nanoplot_filtered.log",
#     benchmark:
#         "{wdir}/{sample}/benchmarks/{genome}.nanoplot_after_filtering.tsv"
#     shell:
#         """
#         NanoPlot --fastq {input.fastq} -t {resources.cpus_per_task} --tsv_stats --outdir {wildcards.wdir}/{wildcards.sample}/nanoplot_filtered/ --prefix '{wildcards.genome}_' --N50 --no_static --verbose --title {wildcards.sample} &> {log}
#         """




