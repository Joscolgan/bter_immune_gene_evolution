#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan)       	Program: run_comp_calc.sh
##
## Date: 28/06/22
##
## Introduction:
## This script takes two input files:
## - A reference genome assembly in FASTA format (provided as an argument)
## - A folder populated by individual BED files corresponding to putative deleted sites.
## The script outputs a tab-delimited file containing each deleted site and its base composition
## as calculated by seqtk 'comp'.
##
#############################################################################################

## Read in input as an argument from the command line:
file=$1

## Index bam file:
bwa index "$file"

## For each BED file (i.e., each putative deletion), perform the following:
## - Create an individual FASTA file for each BED file; 
## - Calculate the proportion of ambiguous (N) bases in the FASTA sequence of each bed.
for name in deletion_bed_files/*bed
do
echo "$name"
new_name="$(echo "$name" | cut -d '.' -f 1 )"
fastaFromBed -fi data/Bombus_terrestris.Bter_1.0.dna.toplevel.fa -bed "$name" -fo "$new_name".fasta
seqtk comp "$new_name".fasta >> deletion_bed_files/combined_comp_calcs.txt;
done
