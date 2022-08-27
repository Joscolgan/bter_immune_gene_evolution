#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

## Read in input FASTA file from the command line:
file=$1

## Create an index using bwa:
bwa index "$file"

