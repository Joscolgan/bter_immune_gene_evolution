#!/usr/bin/env Rscript
##########################################################################################
##
## Author: Joe Colgan (joscolgan)                    Program: nucleotide_div_calculator.R
##
## Date: 20-04-2022
##
## Purpose:
##  - Read in text file containing information of VCF file for each chromosome, as well
##    as start base position and end base position for each corresponding chromosome.
##  - Calculate diversity and neutrality statistics for each chromosome.
##  - Output a table containinig -
##      a) Chromosome
##      b) start position of genomic window
##      c) end position of genomic window
##      d) nucleotide diversity for user-defined genomic window
##      e) number of segregating sites for user-defined genomic window
##      f) Tajima's D for use-defined genomic window
##      g) midpoint position of genomic window
##
##########################################################################################

## Take arguments from the command line:
args = commandArgs(trailingOnly=TRUE)

## Assign 'chromosome_info.txt' to input:
input <- args[1]
print(input)

## Assign reference gff to input_gff:
input_gff <- args[2]
print(input_gff)

## Assign third argument to output
output <- args[3]
print(output)

## Manually set bin width (e.g. size of sliding window):
bin_width <- 100000
jump_width <- bin_width / 2

# Load libraries; install from scratch if needed
libraries <- c("ggplot2",
"PopGenome",
"stringr")
for (lib in libraries) {
    if (require(package = lib, character.only = TRUE)) {
        print("Successful")
    } else {
        print("Installing")
        source("https://bioconductor.org/biocLite.R")
        library(lib, character.only = TRUE )
    }
}

## Turn off scientific notation:
options(scipen = 999)

## Load data
## Read in the tabix indexed VCF
## The input parameters include:
## 1) Tabix-indexed bgzipped file
## 2) numcols: number of SNPs that should be read in as a chunk
## 3) tid: The ID of the chromosome
## 4) Frompos: From position e.g. 1
## 5) topos: End positon (probably length of chromosome)
## An additional parameter is to include a gff

## Scan in information about chromosomes:
chrom_data <- read.table(input, header = FALSE)
chrom_data <- as.matrix(chrom_data)

## Divide each chromosome into sliding windows (window size: 10kb, jump size:5kb):
count <- 0
nucdiv_combined_df <- data.frame()

print(chrom_data)

## For loop, cycle through each chromosomes and calculate nucleotide diverity within
## genomic windows of user-defined size:
for (item in 1:nrow(chrom_data)){
        ## Extract name of chromosome:
        chrom <- gsub(".recode.vcf.bgz",
                      "",
                      chrom_data[item,][1])
        chrom <- gsub("variant_files/",
                      "",
                      chrom)
        ## Read in VCF file:
        chrom_data_item <- readVCF(chrom_data[item,][1],
                                   bin_width,
                                   chrom,
                                   chrom_data[item,][2],
                                   chrom_data[item,][3],
                                   gffpath = input_gff,
                                   approx = FALSE,
                                   include.unknown = TRUE)
        ## Create sliding windows:
        chrom_data_item      <- sliding.window.transform(chrom_data_item, bin_width, jump_width, type=2, whole.data = TRUE)
        ## Extract region names, which contains window coordinates:
        region_names         <- chrom_data_item@region.names
        print(region_names)
        print(str_split_fixed(region_names, " ", 4)[, 1-3])
        region_names_df      <- str_split_fixed(region_names, " ", 4)[, 1-3]
        region_names_df       <- region_names_df[, 1:2]
        print(head(region_names_df))
        ## Calculation of neutrality statistics:
        chrom_data_item      <- neutrality.stats(chrom_data_item)
        ## Calculation of diversity statistics:
        chrom_data_item      <- diversity.stats(chrom_data_item, pi = TRUE)
        ## Calculation of nucleotide diversity:
        nucdiv_item          <- chrom_data_item@nuc.diversity.within
        nucdiv_item          <- nucdiv_item/bin_width
        nucdiv_item.df       <- as.data.frame(cbind(nucdiv_item, region_names_df))
        print(head(nucdiv_item.df))
        nucdiv_item.df$chrom <- chrom
        nucdiv_item.df$tajima_d <- chrom_data_item@Tajima.D
        ## Subset information on the number of segregating sites:
        nucdiv_item.df$seg_sites <- chrom_data_item@n.segregating.sites
        ## Combine nucleotide diversity statistics across windows:
        nucdiv_combined_df <- rbind(nucdiv_combined_df, nucdiv_item.df)
}

print(str(nucdiv_combined_df))

## Rearrange order:
nucdiv_combined_df <- nucdiv_combined_df[c(4,2,3,1,5,6)]

## Rename columns:
colnames(nucdiv_combined_df) <- c("nuc_diversity",
                                  "start",
                                  "end",
                                  "chrom",
                                  "tajima_d",
                                  "seg_sites")

print(head(nucdiv_combined_df))

## Calcule median nucletide diversity
nucdiv_median <- median(as.numeric(as.character(nucdiv_combined_df$nuc_diversity)))

## Export table:
write.table(nucdiv_combined_df,
            output,
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE,
            sep = "\t")
