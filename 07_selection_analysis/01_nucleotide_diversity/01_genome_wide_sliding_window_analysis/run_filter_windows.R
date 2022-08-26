#!/usr/bin/env Rscript
##########################################################################################
##
## Author: Joe Colgan                           Program: run_filter_windows.R
##
## Date: 20-02-2022
##
## Purpose:
##  - Read in text file containing information of VCF file for each chromosome, as well
##    as start base position and end base position for each corresponding chromosome.
##  - Calculate diversity and neutrality statistics for each chromosome.
##  - Output a table containinig -
##	a) Chromosome
##	b) start position of genomic window
##	c) end position of genomic window
##	d) nucleotide diversity for user-defined genomic window
##	e) number of segregating sites for user-defined genomic window
##	f) Tajima's D for use-defined genomic window
##	g) midpoint position of genomic window
##
##########################################################################################

## Take arguments from the command line:
args = commandArgs(trailingOnly=TRUE)

input  <- args[1]
print(input)

output <- args[2]
print(output)

## Read in input:
data <- read.table(file = input,
                   header = FALSE)

## Calculate the percentage of the window with N bases:
data$V14 <- data$V9 / data$V2

## Subset windows where 10% or less of the window consists of N bases:
filtered_data <- subset(x = data,
                        V14 <= 0.1)

## Write file to output:
write.table(x = filtered_data,
            file = output,
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE,
            sep = "\t")
