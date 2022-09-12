#!/usr/bin/env python3
##############################################################################
## 
## Author: Joe Colgan                   Program: calculate_deletion_depth.py
##
## Date: 15/06/2022
##
## Purpose:
## This script takes an alignment BAM file per sample and intersects with multiple 
## BED files containing the genomic positions of individual putative duplication events. 
## The intersection process outputs individual BAM files containing reads aligned
## to each putative duplication (bin). The number of aligned reads per BAM are 
## calculated for each genomic base per individual. The median number of aligned reads
## is calculated and output to a text file for loadiing into R.
##
## This script is a modification of a script previously published by Colgan et
## al. (2022): https://doi.org/10.1093/molbev/msab366
##
##############################################################################
# Import modules
import os.path

from helper_functions import *

# To achieve final output, the Snakefile contains custom-defined rules (see below) that
# outline commands to execute sequentially to take custom-defined input(s) and generate
# final output (as defined in rule all).

# For this specific script, rules are defined as follows:
# rule all:
#   - Defines the expected final output of the Snakefile.
# rule index_genome:
#    - A single fasta file containing genome of interest is provided as input.
#      Input fasta file is indexed by bowtie2.

##############################################################################
# Sample information
##############################################################################
# Bumblebee (Bombus terrestris) males were collected summer 2018
# Each individual was assigned a unique identifier:
# Example: 'MU_01'

##############################################################################
# Prior to use
##############################################################################
# To run calculate_depth.py:
#  1. Download and install the following software:
#   samtools
#   bedtools
#
#  2. Ensure helper_functions.py is within the same directory of the Snakefile.
#
#  3. Assign global variables for use within specific rules
#     Please see section below for further information on variable to be assigned
#
#  4. Assignment of wildcards to be used within the rules
#
#  5. Define all the input and outputs of each rule with respect to the above defined wildcards
#
#  6. Each rule willl take assigned input, global variables and generate defined outputs.
#
#  7. Input data should be formatted in the context of defined wildcards: {samples}_RG_updated.bam
#      For example: MU_01_RG_updated.bam
#
# 8. Make a text.file called 'sample_list.txt' and put in same directory as Snakfile.
#     Populate 'sample_list.txt' with names of samples to be analysed.
#       For example:
#       MU_01
#       MU_02
#       MU_03
#
##############################################################################
# Assign global variables for use in rules (see below)
##############################################################################
# Assign name for list of paths for files containing putative duplications
DUPLICATION_LIST = "CNV_deletions.txt"

##############################################################################
# Assignment of wildcards to be used within rules
##############################################################################
# Open file and read in contents - one sample per line
with open('samples_list.txt') as samples:
    content = samples.readlines()
    SAMPLES = [samples.rstrip('\n') for samples in content]
    print(SAMPLES)

# Read in the content of the duplication list - one duplication per line
with open('CNV_deletions.txt') as duplications:
    dup_content = duplications.readlines()
    DUPLICATIONS = [duplications.rstrip('\n') for duplications in dup_content]
    print(DUPLICATIONS)

##############################################################################
# Specify all input/output files in terms of sample wildcards
##############################################################################
# Assign path for input alignment BAM files.
ALIGNED_DATA            = "results/{samples}_RG_updated.bam"

# Assign path for BED files containing CNV information.
DUPLICATED_DATA         = "filtered_deleted_beds/{duplications}.bed"

# Output intersected BAM files here.
INTERSECTED_DATA        = "raw_deletion_temp/01_intersect_bam/{duplications}.{samples}.bam"

# Output depth counts here.
DEPTH_DATA              = "raw_deletion_temp/02_depth_bam/{duplications}.{samples}.depth.txt"

# Output mean counts per sample here.
CALCULATED_MEDIAN_DATA  = "raw_deletion_temp/03_median_counts/{duplications}.{samples}.median_depth.txt"

# Output combined mean data here.
COMBINED_MEDIAN_DATA    = "raw_deletion_temp/04_combined_counts/{duplications}.combined.median_depth.txt"

# Output combined CNV mean data here.
COMBINED_CNV_DATA       = "raw_deletion_temp/05_combined_CNV_means/combined.median_depth.all_CNVs.txt"

# Output final data file here.
FINAL_DATA              = "raw_deletion_temp/06_final/final_combined.median_depth.all_dups.txt"

##############################################################################
# Define binaries in context of path relative to Snakefile
##############################################################################
# binaries
# Align lines of code using 'Assign Align' using cmd+shift+p and selecting 'align'
# Create dictionaries for directories and  tools'
dirs  = {}
dirs['project']        = os.path.abspath('./')
dirs['src']            = os.path.join(dirs['project'], 'bin')

# Create an empty dictionary called 'tools'
tools = {}
tools['intersect']     = os.path.join(dirs['src'], 'intersectBed')
tools['samtools']      = os.path.join(dirs['src'], 'samtools')

##############################################################################
#
# Specify rules with commands to be executed
#
##############################################################################
# First rule is list the final output
rule all:
    input: expand(FINAL_DATA)

# Each rule will specify an intermediate step
# Intersect BAM file using BED files populated with CNV genomic co-ordinates
rule intersect_bam:
    input:  ALIGNED_DATA, DUPLICATED_DATA
    output: INTERSECTED_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[intersect]} -abam {input[0]} -b {input[1]} > {output} && [[ -s {output} ]]")

# Count depth of coverage for reads aligned within intersected BAM files.
rule count_depth:
    input:  INTERSECTED_DATA
    output: DEPTH_DATA
    run:
        check_files_arent_empty(input)
        shell("{tools[samtools]} depth {input} > {output} && [[ -s {output} ]]")

# Calculate the median of aligned reads for each putative CNV.
rule calculate_median:
    input:  DEPTH_DATA
    output: CALCULATED_MEDIAN_DATA
    run:
        check_files_arent_empty(input)
        shell("Rscript median_counter.R {input} {output} && [[ -s {output} ]]")

# Paste the means for each sample per row
rule paste_means:
    input:  expand("raw_deletion_temp/03_median_counts/{{duplications}}.{samples}.median_depth.txt", samples = SAMPLES)
    output: COMBINED_MEDIAN_DATA
    run:
        check_files_arent_empty(input)
        shell("paste {input} > {output} && [[ -s {output} ]]")

# Combine the CNV mean rows
rule combine_CNVs:
    input:  expand("raw_deletion_temp/04_combined_counts/{duplications}.combined.median_depth.txt", duplications = DUPLICATIONS)
    output: COMBINED_CNV_DATA
    run:
        check_files_arent_empty(input)
        shell("cat {input} > {output} && [[ -s {output} ]]")

# Combine counts with mean row name
rule add_CNV_names:
    input:  COMBINED_CNV_DATA
    output: FINAL_DATA
    run:
        check_files_arent_empty(input)
        shell("paste {DUPLICATION_LIST} {input} > {output} && [[ -s {output} ]]")
