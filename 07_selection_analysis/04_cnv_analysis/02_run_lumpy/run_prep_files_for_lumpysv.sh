#!/bin/sh

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

################################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan).                    Program: run_lumpysv.sh
##
## Date: 22/05/22
## 
## Purpose:
## This script takes for each sample:
## - An alignment file in bam format sorted by position.   
## These files are used by samtools and lumpysv to identify and produce subsetted alignment (bam)
## files containing discordant and split reads, Each bam produced is then sorted by samtools.
##
################################################################################################

for file in results_dtol/*bam;
do 
sample="$(echo $file |cut -d '/' -f 2 | cut -d '.' -f 1)";  
# Extract the discordant paired-end alignments.
samtools view -b -F 1294 results_dtol/"$sample".bam > results_dtol/"$sample".discordants.unsorted.bam

# Extract the split-read alignments
samtools view -h results_dtol/"$sample".bam \
    | src/lumpy-sv-0.3.1/scripts/extractSplitReads_BwaMem -i stdin \
    | samtools view -Sb - \
    > results_dtol/"$sample".splitters.unsorted.bam

# Sort all alignments
samtools sort -@20 -o results_dtol/"$sample".discordants results_dtol/"$sample".discordants.unsorted.bam
samtools sort -@20 -o results_dtol/"$sample".splitters results_dtol/"$sample".splitters.unsorted.bam
samtools sort -@20 -o results_dtol/"$sample".sorted results_dtol/"$sample".bam
done

