#!/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan)       	Program: run_deletions_snakemake.sh
##
## Date: 28/06/22
##
## Introduction:
## The purpose of this script is to call a snakefile to estimate read depth for putative
## deleted sites within the genomes of wild-sampled bumblebee males.
##
#############################################################################################

## Load conda:
module load conda/2

## Activate conda environment:
source activate conda_env

## Load R:
module load r

## Run snakefile:
snakemake -s calculate_deletion_depth.py -p -j 15 --rerun-incomplete
