#!/bin/sh
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

#############################################################################################
## Author: Sarah Larragy, Joe Colgan (joscolgan)         	Program: run_fastp.sh
##
## Date: 30/01/22
##
## Program:
## This script takes FastQC files as input and performs filtering using fastp.
## The output of the script is filtered reads that can be used for genome alignment
#############################################################################################

for file in input/*1.fq.gz;
do
reverse="$(echo $file |cut -d '_' -f 1-5)";
output="$(echo $reverse |cut -d '/' -f 2)" ;
echo $output;
fastp -i $file \
-I "$reverse"_2.fq.gz \
-o results/"$output"_1.fq.gz \
-O results/"$output"_2.fq.gz \
--adapter_sequence AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCT \
--adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG \
--detect_adapter_for_pe \
-q 20 \
-u 40 \
--n_base_limit 15 \
-l 50 \
-p \
--thread 20;
done
