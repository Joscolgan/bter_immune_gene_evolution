#!/bin/sh

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan)       	Program: run_bedtools_depth.sh
##
## Date: 28/06/22
##
## Introduction:
## The purpose of this script is to calculate read depth for putative copy number variation.
## The script takes a folder containing individual BED files corresponding to a CNV, intersects
## the coordinates of each putative CNV with an alignment (bam) file for each individual sample
## and calculates the number of reads covering each base of the putative CNV.
## The script outputs a tab-delimited file for each CNV containing the read depth for each base
## within a CNV.
##
#############################################################################################

for name in bed_files/B*bed;
do
bam_file="$(echo "$name" | cut -d '.' -f 1 )";
echo "$bam_file";
mkdir "$bam_file";
for bam in bed_files/bam_files/MU_*sorted;
do
echo "$bam";
new_bam="$(echo "$bam" | cut -d '/' -f 3 | cut -d '.' -f 1)";
intersectBed -a "$bam" -b "$name" | samtools depth - > "$bam_file"/"$new_bam"_intersect.text;
echo "Done";
done;
done
