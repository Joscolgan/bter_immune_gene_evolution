#!/bin/sh

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

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

