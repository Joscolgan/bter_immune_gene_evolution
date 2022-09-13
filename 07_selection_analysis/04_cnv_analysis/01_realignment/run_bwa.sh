#!/bin/sh

#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

################################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan).                    Program: run_bwa.sh
##
## Date: 20/05/22
## 
## Purpose:
## This script takes compressed fastq files (two per sample; pairs) and aligns against an 
## indexed reference genome assembly using BWA. The script implements a for loop to align 
## each sample in turn and output an alignment file in bam format for each sample.
##
################################################################################################

for file in input/MU_*1.fq.gz;
do
reverse="$(echo $file |cut -d '_' -f 1-5)";
output="$(echo $reverse |cut -d '/' -f 2)" ;
echo $output;
bwa mem -t 40 -R "@RG\tID:$output\tSM:$output\tLB:lib" \
    data/Bombus_terrestris-GCA_000214255.1-unmasked.fa "$file" "$reverse"_2.fq.gz \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > results_dtol/"$output"_dtol_alignment.bam
done
