#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

file=$1

bwa index "$file"

#for file in input/*1.fq.gz;
#do
#reverse="$(echo $file |cut -d '_' -f 1-5)";
#output="$(echo $reverse |cut -d '/' -f 2)" ;
#echo $output;
#bwa mem -t 40 -R "@RG\tID:id\tSM:$output\tLB:lib" data/Bombus_terrestris.Bter_1.0.dna.toplevel.fa "$file" "$reverse"_2.fq.gz \
#    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
#    | samtools view -S -b - \
#    > results/"$output".bam
#done