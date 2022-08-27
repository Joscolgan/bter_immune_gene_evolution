#!/bin/sh

#SBATCH --time=10:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

for file in input/MU_2*1.fq.gz;
do
reverse="$(echo $file |cut -d '_' -f 1-5)";
output="$(echo $reverse |cut -d '/' -f 2)" ;
echo $output;
bwa mem -t 40 -R "@RG\tID:$output\tSM:$output\tLB:lib" data/Bombus_terrestris-GCA_000214255.1-unmasked.fa "$file" "$reverse"_2.fq.gz \
    | samblaster --excludeDups --addMateTags --maxSplitCount 2 --minNonOverlap 20 \
    | samtools view -S -b - \
    > results_dtol/"$output"_dtol_alignment.bam
done
