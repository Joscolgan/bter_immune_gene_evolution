#!/bin/sh

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

for name in shuffle_analysis/B01*bed;
do
bam_file="$(echo "$name" | cut -d '.' -f 1 )";
echo "$bam_file";
mkdir "$bam_file";
for bam in shuffle_analysis/bam_files/MU_*sorted;
do
echo "$bam";
new_bam="$(echo "$bam" | cut -d '/' -f 3 | cut -d '.' -f 1)";
intersectBed -a "$bam" -b "$name" | samtools depth - > "$bam_file"/"$new_bam"_intersect.text;
echo "Done";
done;
done