#!/bin/sh

#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

module load conda/2
source activate conda_env

module load parallel

ref="data/Bombus_terrestris.Bter_1.0.dna.toplevel.fa"

src/freebayes/scripts/freebayes-parallel \
region_list.txt 40 \
-f "$ref" \
--bam-list irish_bam_list.txt \
--ploidy 2 \
--report-genotype-likelihood-max \
--use-mapping-quality \
--genotype-qualities \
--use-best-n-alleles 4 \
--haplotype-length 0 \
--min-base-quality 3 \
--min-mapping-quality 1 \
--min-alternate-frac 0.25 \
--min-coverage 1 \
--use-reference-allele > results/irish_samples.vcf
