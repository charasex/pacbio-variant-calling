rule mosdepth_summary:
    input:
        bam="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam",
        bai="{wdir}/{sample}/bam/{genome}_minimap2_sorted.bam.bai"
    output:
        dist="{wdir}/{sample}/mapping_QC/{genome}_mosdepth.mosdepth.global.dist.txt",
        summary="{wdir}/{sample}/mapping_QC/{genome}_mosdepth.mosdepth.summary.txt",
        coverage_windows="{wdir}/{sample}/mapping_QC/{genome}_mosdepth.regions.bed.gz",
        coverage_windows_idx="{wdir}/{sample}/mapping_QC/{genome}_mosdepth.regions.bed.gz.csi"
    conda:
        "../envs/mosdepth.yaml"
    threads: 8
    log:
        "{wdir}/{sample}/logs/{genome}_mosdepth.log"
    benchmark:
        "{wdir}/{sample}/benchmarks/{genome}_mosdepth_summary.tsv"
    params:
        mosdepth_windows_size=config.get("mosdepth_windows_size", 10000),
        prefix="{wdir}/{sample}/mapping_QC/{genome}_mosdepth"
    shell:
        """
        mosdepth --no-per-base \
            -t {threads} \
            --by {params.mosdepth_windows_size} \
            {params.prefix} \
            {input.bam} 2> {log}
        """