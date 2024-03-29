#!/bin/sh

#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

################################################################################################
##
## Author: Sarah Larragy, Joe Colgan (joscolgan).                    Program: run_lumpysv.sh
##
## Date: 22/05/22
## 
## Purpose:
## This script takes for each sample:
## - An alignment file in bam format sorted by position.  
## - An alignment file in bam format populated with discordant reads.  
## - An alignment file in bam format populated with split reads.   
## These files are used by lumpysv to identify putative copy number variants present in the input
## population. The script outputs a file in variant call format (VCF) providing information on
## the physical genomic locations of putative sites, as well as measures of support, as well
## as genotype calls for each individual.
##
################################################################################################

lumpyexpress -B \
results_dtol/MU_01_FDPL190630198-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_03_FDPL190630200-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_06_FDPL190630203-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_07_FDPL190630204-1a_HKYWTDSXX_L24_dtol_alignment.sorted,\
results_dtol/MU_08_FDPL190630205-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_09_FDPL190630206-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_10_FDPL190630207-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_12_FDPL190630209-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_17_FDPL190630214-1a_HHYHCDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_18_FDPL190630215-1a_HHYHCDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_19_FDPL190630216-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_20_FDPL190630217-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_21_FDPL190630218-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_22_FDPL190630219-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_23_FDPL190630220-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_25_FDPL190630222-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_26_FDPL190630223-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_27_FDPL190630224-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_28_FDPL190630225-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_29_FDPL190630226-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_32_FDPL190630229-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_33_FDPL190630230-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_35_FDPL190630231-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_36_FDPL190630232-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_38_FDPL190630233-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_39_FDPL190630234-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_40_FDPL190630235-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_41_FDPL190630236-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_42_FDPL190630237-1a_HJL3TDSXX_L2_dtol_alignment.sorted,\
results_dtol/MU_43_FDPL190630238-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_44_FDPL190630239-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_45_FDPL190630240-1a_HJL3TDSXX_L1_dtol_alignment.sorted,\
results_dtol/MU_46_FDPL190630241-1a_HL272DSXX_L1_dtol_alignment.sorted \
-D \
results_dtol/MU_01_FDPL190630198-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_03_FDPL190630200-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_06_FDPL190630203-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_07_FDPL190630204-1a_HKYWTDSXX_L24_dtol_alignment.discordants,\
results_dtol/MU_08_FDPL190630205-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_09_FDPL190630206-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_10_FDPL190630207-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_12_FDPL190630209-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_17_FDPL190630214-1a_HHYHCDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_18_FDPL190630215-1a_HHYHCDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_19_FDPL190630216-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_20_FDPL190630217-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_21_FDPL190630218-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_22_FDPL190630219-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_23_FDPL190630220-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_25_FDPL190630222-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_26_FDPL190630223-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_27_FDPL190630224-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_28_FDPL190630225-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_29_FDPL190630226-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_32_FDPL190630229-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_33_FDPL190630230-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_35_FDPL190630231-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_36_FDPL190630232-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_38_FDPL190630233-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_39_FDPL190630234-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_40_FDPL190630235-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_41_FDPL190630236-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_42_FDPL190630237-1a_HJL3TDSXX_L2_dtol_alignment.discordants,\
results_dtol/MU_43_FDPL190630238-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_44_FDPL190630239-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_45_FDPL190630240-1a_HJL3TDSXX_L1_dtol_alignment.discordants,\
results_dtol/MU_46_FDPL190630241-1a_HL272DSXX_L1_dtol_alignment.discordants \
-S \
results_dtol/MU_01_FDPL190630198-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_03_FDPL190630200-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_06_FDPL190630203-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_07_FDPL190630204-1a_HKYWTDSXX_L24_dtol_alignment.splitters,\
results_dtol/MU_08_FDPL190630205-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_09_FDPL190630206-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_10_FDPL190630207-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_12_FDPL190630209-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_17_FDPL190630214-1a_HHYHCDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_18_FDPL190630215-1a_HHYHCDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_19_FDPL190630216-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_20_FDPL190630217-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_21_FDPL190630218-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_22_FDPL190630219-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_23_FDPL190630220-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_25_FDPL190630222-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_26_FDPL190630223-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_27_FDPL190630224-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_28_FDPL190630225-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_29_FDPL190630226-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_32_FDPL190630229-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_33_FDPL190630230-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_35_FDPL190630231-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_36_FDPL190630232-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_38_FDPL190630233-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_39_FDPL190630234-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_40_FDPL190630235-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_41_FDPL190630236-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_42_FDPL190630237-1a_HJL3TDSXX_L2_dtol_alignment.splitters,\
results_dtol/MU_43_FDPL190630238-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_44_FDPL190630239-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_45_FDPL190630240-1a_HJL3TDSXX_L1_dtol_alignment.splitters,\
results_dtol/MU_46_FDPL190630241-1a_HL272DSXX_L1_dtol_alignment.splitters \
-o variant_calls_raw_dtol.vcf \
-P
