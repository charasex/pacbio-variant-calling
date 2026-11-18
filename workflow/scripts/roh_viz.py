import os
import sys
import traceback
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Rediriger stderr vers le log Snakemake si défini
if snakemake.log:
    sys.stderr = open(snakemake.log[0], "w")

try:
    df = pd.read_csv(
        snakemake.input.win_stats,
        sep="\t",
        names=["chrom", "start", "end", "dp", "snps", "het_rate"]
    )

    if os.path.exists(snakemake.input.roh) and os.path.getsize(snakemake.input.roh) > 0:
        roh_df = pd.read_csv(
            snakemake.input.roh,
            sep="\t",
            names=["chrom", "start", "end", "length"]
        )
    else:
        roh_df = pd.DataFrame(columns=["chrom", "start", "end", "length"])

    # Sélection des principaux contigs (plafonné aux 25 plus grands pour éviter de saturer Plotly)
    chrom_sizes = df.groupby("chrom")["end"].max().sort_values(ascending=False)
    main_chroms = chrom_sizes[chrom_sizes >= snakemake.params.min_contig_len].head(25).index.tolist()
    if not main_chroms:
        main_chroms = chrom_sizes.head(10).index.tolist()

    df = df[df["chrom"].isin(main_chroms)]

    n_rows = max(1, len(main_chroms))
    # Sécurité Plotly : vertical_spacing doit être strictement < 1 / (n_rows - 1)
    v_space = min(0.03, 0.5 / max(1, n_rows - 1)) if n_rows > 1 else 0.02

    fig = make_subplots(
        rows=n_rows,
        cols=1,
        subplot_titles=main_chroms if main_chroms else ["Aucune donnée"],
        vertical_spacing=v_space,
        specs=[[{"secondary_y": True}] for _ in range(n_rows)]
    )

    for idx, chrom in enumerate(main_chroms, start=1):
        sub = df[df["chrom"] == chrom]
        mid_pos = (sub["start"] + sub["end"]) / 2

        fig.add_trace(
            go.Scatter(
                x=mid_pos,
                y=sub["het_rate"] * 1000,
                mode="lines",
                name="SNPs 0/1 / kb",
                line=dict(width=1.2, color="#1f77b4"),
                showlegend=(idx == 1)
            ),
            row=idx, col=1, secondary_y=False
        )

        fig.add_trace(
            go.Scatter(
                x=mid_pos,
                y=sub["dp"],
                mode="lines",
                name="Profondeur (DP)",
                line=dict(width=0.8, color="rgba(120,120,120,0.45)"),
                showlegend=(idx == 1)
            ),
            row=idx, col=1, secondary_y=True
        )

        for _, r in roh_df[roh_df["chrom"] == chrom].iterrows():
            fig.add_vrect(
                x0=r["start"],
                x1=r["end"],
                fillcolor="red",
                opacity=0.25,
                line_width=0,
                row=idx, col=1
            )

        fig.update_yaxes(title_text="SNPs/kb", row=idx, col=1, secondary_y=False)
        fig.update_yaxes(title_text="DP", row=idx, col=1, secondary_y=True, showgrid=False)

    fig.update_layout(
        height=240 * n_rows,
        title=f"Hétérozygotie (SNPs/kb), Couverture et ROH — {snakemake.wildcards.sample} ({snakemake.wildcards.primary})",
        template="plotly_white"
    )

    fig.write_html(snakemake.output.roh_plot)

except Exception:
    traceback.print_exc()
    sys.exit(1)