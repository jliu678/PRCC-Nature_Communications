#!/bin/bash

#SBATCH --partition=bigmem
#SBATCH --job-name=snakemake_single
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --output=/data/wanglab_mgberis/OUT.%x.%j
#SBATCH --error=/data/wanglab_mgberis/ERR.%x.%j
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

source /PHShome/${USER}/.bashrc
conda activate snakemake

#   Load a suitable module for snakemake
module load STAR

code_dir=/data/wanglab_mgberis/siyi2023summer_snakemake_10x
#   cd
mkdir /data/wanglab_mgberis/tmp
cd /data/wanglab_mgberis/tmp

#   run snakemake
#snakemake --cores 12 --configfile ${code_dir}/config_10x5p.yaml --snakefile ${code_dir}/Snakefile_10x5prime
snakemake --cores 24 --configfile ${code_dir}/config_10x3v2or3.yaml --snakefile ${code_dir}/Snakefile_10x3primev2or3