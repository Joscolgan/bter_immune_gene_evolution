#!/bin/sh

#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

################################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan).                    Program: make_bwa_index.sh
##
## Date: 20/05/22
## 
## Purpose:
## This script takes a reference genome assembly in FASTA format from the command line and
## with BWA creates an index for genome alignment.
##
################################################################################################

## Read in input FASTA file from the command line:
file=$1

## Create an index using bwa:
bwa index "$file"
