#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
## Author: Sarah Larragy, Joe Colgan (joscolgan)                 Program: run_bowtie2.sh
##
## Date: 28/06/22
##
## Introduction:
## The purpose of this script is to take fastp-filtered FASTQ sequences and align
## against a reference genome assembly using the short read aligner Bowtie2.
## The script generates a compressed alignment file (BAM format) for each samples.
## The output files are sorted using samtools.
## 
#############################################################################################


for file in input/*1.fq.gz;
do
reverse="$(echo "$file" |cut -d '_' -f 1-5)_2.fq.gz";
output="$(echo "$reverse" |cut -d '/' -f 2 | cut -d '_' -f 1-5)";
echo "$output";
bowtie2 --rg-id "$output" \
--rg SM:"$output" \
--rg LB:library1 \
--rg PL:ILLUMINA \
--rg DS:Novaseq6000 \
--local -p 40 \
-X 1000 -x data/genome/bter_1.0 \
-1 <(zcat "$file") \
-2 <(zcat "$reverse") | \
samtools view -bS | samtools sort -o results/"$output".bam -;
done
