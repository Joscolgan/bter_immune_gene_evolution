#!/bin/sh

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
## 
## Author: Sarah Larragy, Joe Colgan (joscolgan)         Program: convert_vcf_to_admixture.sh
##
## Date: 15/04/2022
##
## Introduction:
## The purpose of this script is to take freebayes-called and vcftools-filtered variant
## calls and generated input files for analysis with ADMIXTURE.
## The output of the script is a vcf, which is used as input into ADMIXTURE.
## 
#############################################################################################

## take input and output arguments from the command line:
vcf=$1
positions=$2
output=$3

## Extract positions for running admixture:
vcftools --vcf "$vcf" \
--positions "$positions" \
--recode \
--out "$output". 
