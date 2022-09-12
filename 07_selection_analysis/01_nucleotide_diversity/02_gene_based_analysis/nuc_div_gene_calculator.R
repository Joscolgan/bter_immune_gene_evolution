#!/usr/bin/env Rscript
##########################################################################################
##
## Author: Joe Colgan (joscolgan)                      Program: nuc_div_gene_calculator.R
##
## Date: 30-04-2020
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

input <- "chromosome_info.txt"
bin_width <- 10000
input_gff <- "./gff/Bombus_terrestris.Bter_1.0.51.gff3"
output <- "../results/theta_estimates_per_gene.txt"

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

## Turn off scientific notations:
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

## Define jumpwidth:
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
        ## Adds the results to the appropriate slots
        chrom_data_item <- splitting.data(chrom_data_item,
                                          subsites = "gene")
        chrom_data_item <- diversity.stats(chrom_data_item)
        print(chrom_data_item@region.names)
        # Extract region names, which contains window coordinates:
        region_names         <- chrom_data_item@region.names
        print(region_names)
        print(str_split_fixed(region_names, " ", 4)[, 1-3])
        region_names_df      <- str_split_fixed(region_names, " ", 4)[, 1-3]
        region_names_df      <- region_names_df[, 1:2]
        gene_start           <- as.numeric(region_names_df[, 1])
        gene_end             <- as.numeric(region_names_df[, 2])
        gene_length          <- gene_end - gene_start
        print(head(region_names_df))
        print(class(region_names_df))
        ## Calculation of neutrality statistics:
        chrom_data_item      <- neutrality.stats(chrom_data_item)
        ## Calculation of diversity statistics:
        chrom_data_item      <- diversity.stats(chrom_data_item,
                                                pi = TRUE)
        ## Calculation of nucleotide diversity:
        nucdiv_item          <- chrom_data_item@nuc.diversity.within
        nucdiv_item.df       <- as.data.frame(cbind(nucdiv_item,
                                                    region_names_df))
        colnames(nucdiv_item.df) <- c("nuc_diversity",
                                      "start",
                                      "end")
        print(head(nucdiv_item.df))
        nucdiv_item.df$chrom <- chrom
        nucdiv_item.df$tajima_d <- chrom_data_item@Tajima.D
        ## Add gene length:
        nucdiv_item.df$gene_length <- gene_length
        ## Adjust nucleotide diversity by gene length:
        nucdiv_item.df$test <- as.numeric(nucdiv_item.df$nuc_diversity) / nucdiv_item.df$gene_length
        ## Combine nucleotide diversity statistics across windows:
        nucdiv_combined_df <- rbind(nucdiv_combined_df,
                                    nucdiv_item.df)
}

print("Complete")
dim(nucdiv_combined_df)

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
