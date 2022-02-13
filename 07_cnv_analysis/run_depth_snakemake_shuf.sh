#!/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

module load conda/2

source activate conda_env


snakemake -s calculate_duplication_depth_shuf.py -p -j 15
