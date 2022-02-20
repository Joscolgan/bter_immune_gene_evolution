#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

module load conda/2

source activate conda_env

module load r

snakemake -s calculate_depth.py -p -j 10 --rerun-incomplete
