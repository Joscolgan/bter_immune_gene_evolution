#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan)       	Program: run_duplications_snakemake.sh
##
## Date: 28-06-22
##
## Introduction:
## The purpose of this script is to call a snakefile to estimate read depth for putative
## duplicated sites within the genomes of wild-sampled bumblebee males.
##
#############################################################################################

## Load conda:
module load conda/2

## Activate conda environment:
source activate conda_env

## Load R:
module load r

## Run snakemake:
snakemake -s calculate_duplication_depth.py -p -j 15 --rerun-incomplete
