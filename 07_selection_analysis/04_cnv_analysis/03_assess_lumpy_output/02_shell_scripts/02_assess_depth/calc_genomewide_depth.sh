#!/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan)       	Program: calc_genomewide_depth.sh
##
## Date: 28/06/22
##
## Introduction:
## The purpose of this script is to take sorted alignment (bam) files as input and calculate
## depth (base coverage) for each site in the genomes of wild-sampled bumblebees. 
## For each sample, the script outputs a tab-delimited text file containing the physical
## location of the base (Chromosome, base position) and the corresponding depth (i.e., how
## many reads cover a particular site in the genome).
##
#############################################################################################

for name in results/*sorted
do
echo "$name";
new_name="$(echo "$name" | cut -d '.' -f 1 )"
samtools depth "$name" > "$new_name".depth.txt
done
