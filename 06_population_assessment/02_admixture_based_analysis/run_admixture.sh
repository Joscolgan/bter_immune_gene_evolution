#!/bin/sh

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
## 
## Author: Sarah Larragy, Joe Colgan (joscolgan)             Program: run_admixture.sh
##
## Date: 16/04/22
##
## Introduction:
## The purpose of this script is to take a plink.bed file generated from a filtered VCF
## containing positions of polymorphic sites in a natural bumblebee population. This script
## run the software tool ADMIXTURE to determine co-ancestry shared amongst samples in the 
## input data file. The software performs cross-validation for each tested population
## and outputs a log file for each level of K tested.
##
#############################################################################################

for K in {1..20};
do
admixture --cv results/test.plink.bed \
-j20 $K | tee results/log${K}.out;
done
