#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

for file in input/*1.fq.gz;
do
reverse="$(echo "$file" |cut -d '_' -f 1-5)_2.fq.gz";
output="$(echo "$reverse" |cut -d '/' -f 2 | cut -d '_' -f 1-5)";
echo "$output";
bowtie2 --rg-id "$output" \
--rg SM:"$output" \
--rg LB:library1 \
--rg PL:ILLUMINA \
--rg DS:HiSeq2500 \
--local -p 40 \
-X 1000 -x data/genome/bterr_1.0 \
-1 <(zcat "$file") \
-2 <(zcat "$reverse") | \
samtools view -bS | samtools sort -o results/"$output".bam -;
done
