# rule roh:
#     """
#     Détecte les ROH sur primary_haplotype à partir de la chute de densité des SNPs 0/1
#     en vérifiant que la couverture mosdepth est suffisante.
#     """
#     input:
#         vcf_clean_snp="{wdir}/{sample}/mpileup/{primary}_bcftools_clean_snp.vcf.gz",
#         mosdepth_bed="{wdir}/{sample}/mapping_QC/{primary}_mosdepth.regions.bed.gz",
#         mosdepth_summary="{wdir}/{sample}/mapping_QC/{primary}_mosdepth.mosdepth.summary.txt"
#     output:
#         win_stats="{wdir}/{sample}/roh/{primary}_het_windows.tsv",
#         roh="{wdir}/{sample}/roh/{primary}_roh_results.bed"
#     params:
#         min_dp=config["DP"],
#         max_het_rate=config["max_het_rate"],   
#         merge_dist=config["merge_dist"],      
#         min_roh_len=config["min_roh_len"]     
#     conda:
#         "../envs/bedtools.yaml"
#     log:
#         "{wdir}/{sample}/logs/{primary}_roh.log"
#     shell:
#         """
#         (
#         MEAN_DP=$(awk '$1=="total" {{print $4}}' {input.mosdepth_summary})

#         # 1. Compter les SNPs 0/1 par fenêtre mosdepth
#         bedtools intersect -a {input.mosdepth_bed} -b {input.vcf_clean_snp} -c | \
#         awk 'BEGIN{{OFS="\\t"}} {{len=$3-$2; rate=(len>0 ? $5/len : 0); print $1,$2,$3,$4,$5,rate}}' > {output.win_stats}

#         # 2. Garder les fenêtres diploïdes (0.6x <= DP <= 1.8x) à faible hétérozygotie et fusionner
#         awk -v mean_dp="$MEAN_DP" -v max_het={params.max_het_rate} \
#             '$4 >= (0.6 * mean_dp) && $4 <= (1.8 * mean_dp) && $6 <= max_het {{print $1"\\t"$2"\\t"$3}}' {output.win_stats} | \
#         bedtools merge -i - -d {params.merge_dist} | \
#         awk -v min_len={params.min_roh_len} '($3-$2) >= min_len {{print $0"\\t"($3-$2)}}' > {output.roh}
#         ) 2> {log}
#         """

rule bcftools_roh_pass:
    input:
        vcf_all="{wdir}/{sample}/mpileup/{primary}_bcftools_allsites.vcf.gz",
        csi_all="{wdir}/{sample}/mpileup/{primary}_bcftools_allsites.vcf.gz.csi",
        fai="{wdir}/genome/{primary}.fasta.fai"
    output:
        vcf_pass=temp("{wdir}/{sample}/mpileup/{primary}_bcftools_allsites_PASS.vcf.gz"),
        csi_pass=temp("{wdir}/{sample}/mpileup/{primary}_bcftools_allsites_PASS.vcf.gz.csi"),
        roh_full="{wdir}/{sample}/roh/{primary}_bcftools_roh_full.txt.gz",
        roh_rg="{wdir}/{sample}/roh/{primary}_bcftools_roh_RG.txt",
        roh_bed="{wdir}/{sample}/roh/{primary}_bcftools_roh.bed"
    params:
        AFdflt=config["AFdefault"],
        G=config["G"],
        min_chr_len=5000000,  
        min_qual=30
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{primary}_bcftools_roh.log"
    shell:
        """
        (
        CHROMS=$(awk -v min={params.min_chr_len} '$2 >= min {{print $1}}' {input.fai} | paste -sd, -)

        # 1. Fichier PASS restreint aux 31 chromosomes
        bcftools view -f PASS -r "$CHROMS" -Oz -o {output.vcf_pass} {input.vcf_all}
        bcftools index {output.vcf_pass}

        # 2. Sortie complète (ST + RG) compressée pour roh-viz
        bcftools roh --AF-dflt {params.AFdflt} -G {params.G} {output.vcf_pass} | \
            gzip -c > {output.roh_full}

        # 3. Extraction des régions RG et conversion en BED
        zgrep -E "^#|^RG" {output.roh_full} > {output.roh_rg}
        awk -v q={params.min_qual} 'BEGIN{{OFS="\\t"}} $1=="RG" && $8>=q {{print $3, $4-1, $5, $6, $8}}' \
            {output.roh_rg} > {output.roh_bed}
        ) 2> {log}
        """

rule bcftools_roh_viz:
    input:
        roh_full="{wdir}/{sample}/roh/{primary}_bcftools_roh_full.txt.gz",
        vcf_pass="{wdir}/{sample}/mpileup/{primary}_bcftools_allsites_PASS.vcf.gz",
        csi_pass="{wdir}/{sample}/mpileup/{primary}_bcftools_allsites_PASS.vcf.gz.csi",
        fai="{wdir}/genome/{primary}.fasta.fai"
    output:
        html="{wdir}/{sample}/roh/{primary}_roh_viz.html"
    params:
        min_chr_len=5000000
    conda:
        "../envs/bcftools.yaml"
    log:
        "{wdir}/{sample}/logs/{primary}_roh_viz.log"
    shell:
        """
        (
        CHROMS=$(awk -v min={params.min_chr_len} '$2 >= min {{print $1}}' {input.fai} | paste -sd, -)

        roh-viz -i {input.roh_full} -v {input.vcf_pass} -r "$CHROMS" -o {output.html}
        ) 2> {log}
        """


# rule roh_viz:
#     """
#     Visualise la densité de SNPs hétérozygotes, la profondeur et les blocs ROH en HTML interactif.
#     """
#     input:
#         win_stats="{wdir}/{sample}/roh/{primary}_het_windows.tsv",
#         roh="{wdir}/{sample}/roh/{primary}_roh_results.bed"
#     output:
#         roh_plot="{wdir}/{sample}/roh/{primary}_roh_plot.html"
#     params:
#         min_contig_len=config["min_contig_len"]
#     conda:
#         "../envs/plotly.yaml"
#     log:
#         "{wdir}/{sample}/logs/{primary}_roh_viz.log"
#     script:
#         "../scripts/roh_viz.py"