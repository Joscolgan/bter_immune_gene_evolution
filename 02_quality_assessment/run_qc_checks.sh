#!/env/sh
#############################################################################################
## 
## Author: Sarah Larragy, Joe Colgan (github:joscolgan)	       Program: run_qc_checks.sh
##
## Date: 28/01/22
##
## Introduction:
## The purpose of the script is to quality assess FASTQ sequences by:
## - Running FastQC on each sample.
## - Running MultiQC on the FastQC output files to generate a summary.
##
#############################################################################################

## 1. Create project directories:
mkdir fastqc
cd fastqc

## Create subdirectories:
mkdir input
mkdir results
mkdir src

## 2. Locally install software tools:
module load conda/2
module load java
conda create -n conda_env
source activate conda_env
conda install -c bioconda -c conda-forge multiqc

cd src
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip

## Ensure FastQC is executable:
chmod +X fastqc

## Link to local bin:
cd ../bin/
ln -s ../fastqc_project/src/FastQC .
cd ../

## Linked raw, zipped genome files for each sample to /fastqc/input directory:
cd input
ln -s ../../data/X201SC19031072-Z01-F001/raw_data/*/*.gz .

## On DevQ:
srun -p DevQ -N 1 -A mulif005c -t 1:00:00 --pty bash

#Once all files are in the same /input directory, can run fastqc in the same /input directory:
../bin/fastqc -t 10 *gz

## To run MultiQC on all files in current directory:
multiqc .
