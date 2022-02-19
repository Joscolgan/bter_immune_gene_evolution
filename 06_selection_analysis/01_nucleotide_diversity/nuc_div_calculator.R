#!/usr/bin/env Rscript
##########################################################################################
##
## Author: Joe Colgan                           Program: nuc_div_calculator.R
##
## Date: 04-03-2018
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
bin_width <- 100000
input_gff <- "../gff/Bombus_terrestris.Bter_1.0.45.gff"
output <- "bter_hib_all_chrom_100kb.txt"
pop1 <- scan(file = "../populations/irish_samples.txt", as.character())
pop2 <- scan(file = "../populations/non_irish_samples.txt", as.character())

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

options(scipen=999)

#install.packages("WhopGenome")

## Load data
## Read in the tabix indexed VCF
## The input parameters include:
## 1) Tabix-indexed bgzipped file
## 2) numcols: number of SNPs that should be read in as a chunk
## 3) tid: The ID of the chromosome e.g. "NC_015764.1"
## 4) Frompos: From position e.g. 1
## 5) topos: End positon (probably length of chromosome)
## An additional parameter is to include a gff

## Scan in information about chromosomes:
chrom_data <- read.table(input, header = FALSE)
chrom_data <- as.matrix(chrom_data)

## Divide each chromosome into sliding windows (window size: 10kb, jump size:5kb)
##
count <- 0
nucdiv_combined_df <- data.frame()
fst_estimates_df <- data.frame()

## Define jumpwidth:
jump_width <- as.numeric(bin_width)/2
print(jump_width)
print(chrom_data)

## For loop, cycle through each chromosomes and calculate nucleotide diverity within
## genomic windows of user-defined size:
for (item in 1:nrow(chrom_data)){
        ## Extract name of chromosome:
        chrom <- gsub(".recode.vcf.bgz", "", chrom_data[item,][1])
        ## Read in VCF file:
        chrom_data_item <- readVCF(chrom_data[item,][1], bin_width, chrom, chrom_data[item,][2], chrom_data[item,][3], gffpath = input_gff, approx = TRUE,  include.unknown=TRUE)
        ## Create sliding windows:
        chrom_data_item      <- sliding.window.transform(chrom_data_item, bin_width, jump_width, type=2, whole.data = TRUE)
        ## Extract region names, which contains window coordinates:
        region_names         <- chrom_data_item@region.names
        print(region_names)
        print(str_split_fixed(region_names, " ", 4)[, 1-3])
        #region_names_df      <- str_split_fixed(region_names, " ", 4)[,4]
        region_names_df      <- str_split_fixed(region_names, " ", 4)[, 1-3]
        region_names_df       <- region_names_df[, 1:2]
	#region_names_df      <- gsub(" :", "", region_names_df)
        #region_names_df      <- str_split_fixed(region_names_df, " - ", 2)
        print(head(region_names_df))
        ## Calculation of neutrality statistics:
        chrom_data_item      <- neutrality.stats(chrom_data_item)
        ## Calculation of diversity statistics:
        chrom_data_item      <- diversity.stats(chrom_data_item, pi = TRUE)
        ## Calculation of nucleotide diversity:
        nucdiv_item          <- chrom_data_item@nuc.diversity.within
        nucdiv_item          <- nucdiv_item/bin_width
        nucdiv_item.df <- as.data.frame(cbind(nucdiv_item, region_names_df))
        print(head(nucdiv_item.df))
        nucdiv_item.df$chrom <- chrom
        nucdiv_item.df$tajima_d <- chrom_data_item@Tajima.D
        ## Subset information on the number of segregating sites:
        nucdiv_item.df$seg_sites <- chrom_data_item@n.segregating.sites
        ## Combine nucleotide diversity statistics across windows:
        nucdiv_combined_df <- rbind(nucdiv_combined_df, nucdiv_item.df)
        ## Set populations
        chrom_data_item <- set.populations(chrom_data_item, list(pop1, pop2))
        # Diversities and FST (by scaffold)
        chrom_data_item <- F_ST.stats(chrom_data_item) # this does the calculations and 
                       # adds the results to the appropriate slots

# Print FST
fst_estimates <- as.data.frame(get.F_ST(chrom_data_item)) # each line is a scaffold
fst_estimates_df <- rbind(fst_estimates_df, fst_estimates)
chrom_data_item@nucleotide.F_ST
}

## Rename columns:
colnames(nucdiv_combined_df) <- c("nuc_diversity",
                                  "start",
                                  "end",
                                  "chrom",
                                  "tajima_d",
                                  "seg_sites")

print(head(nucdiv_combined_df))

## Rearrange order:
nucdiv_combined_df <- nucdiv_combined_df[c(4,2,3,1,5,6)]

## Calcule median nucletide diversity
nucdiv_median <- median(as.numeric(as.character(nucdiv_combined_df$nuc_diversity)))

## Convert to numeric values:
#nucdiv_combined_df$nuc_diversity <- as.numeric(as.character(nucdiv_combined_df$nuc_diversity))
##nucdiv_combined_df$start         <- as.numeric(as.character(unlist(nucdiv_combined_df$start)))
##nucdiv_combined_df$end           <- as.numeric(as.character(unlist(nucdiv_combined_df$end)))
#nucdiv_combined_df$seg_sites     <- as.numeric(as.character(nucdiv_combined_df$seg_sites))
#nucdiv_combined_df$tajima_d      <- as.numeric(as.character(nucdiv_combined_df$tajima_d))
#nucdiv_combined_df$midpoint      <- round((nucdiv_combined_df$start + nucdiv_combined_df$end)/2)

## Export table:
write.table(nucdiv_combined_df,
output,
row.names = FALSE,
col.names = TRUE,
quote = FALSE,
sep="\t")

## Export table:
write.table(fst_estimates_df,
"./fst_estimates.txt",
row.names = FALSE,
col.names = TRUE,
quote = FALSE,
sep="\t")
