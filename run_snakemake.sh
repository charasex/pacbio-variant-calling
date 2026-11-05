#!/bin/bash
#
#SBATCH --mail-user=ishaac.chninak@inrae.fr
#SBATCH --mail-type=all
#SBATCH --cpus-per-task=10
#SBATCH --mem=20GB
#SBATCH --account=ichninak
#SBATCH --time=20:00:00
#SBATCH --job-name=Charasex
#SBATCH --output=log/slurm-%A.out

species=$1

# module load snakemake/8.9.0
# module load snakemake/8.27.1
# module load snakemake/9.4.0
# . /local/env/envconda.sh
# source activate snakemake
# module load conda
# source activate snakemake


# echo "Snakemake version"
# snakemake --version
# echo "SLURM version"
# srun --version
# echo "Conda version"
# conda --version

# git status

# mkdir --parents log

# echo "Running Snakemake pipeline for species $species..."

#snakemake -s workflow/Snakefile --configfile data-lewontin/config/config.yaml \
#--use-conda --conda-frontend conda --profile ./profiles/slurm \
#--rerun-incomplete --printshellcmds --unlock \
#--config samples="data-lewontin/config/samples_$species.tsv"

#snakemake -s workflow/Snakefile --configfile data-lewontin/config/config.yaml \
#--use-conda --conda-frontend conda --profile ./profiles/slurm \
#--rerun-incomplete --printshellcmds --slurm-status-attempts 100000 \
#--config samples="data-lewontin/config/samples_$species.tsv"

snakemake --use-conda --conda-frontend mamba --printshellcmds --cores 10
