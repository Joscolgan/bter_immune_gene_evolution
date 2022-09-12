#!/bin/sh

#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ


## Load conda:
module load conda/2

## Activate conda environment:
source activate conda_env

## Load R:
module load r

## Run snakemake 
snakemake -s calculate_depth.py -p -j 10 --rerun-incomplete
