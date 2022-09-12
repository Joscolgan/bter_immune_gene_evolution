#!/bin/sh

################################################################################################
##
## Author: Joe Colgan (joscolgan).                          Program: window_cleaner.sh
##
## Date: 20-04-2022
## 
## Purpose:
## A script to take a tab-delimited text file as input:
## - Extract the gene/window coordinates;
## - Count the number of ambigious bases in each window;
## - Remove windows where bases make up at least 10 percent as may be the result of gaps/technical
## artefacts;
## - Run zscore analysis to identify windows with sig. higher or lower genetic diversity compared 
## to rest of genome:
## The script outputs tab-delimited text files of genes of interest found in regions of elevated
## or reduced nucleotide diversity:
##
################################################################################################

## Load R:
module load r

## Take input and output as arguments from the command line:
input=$1
reference=$2
gff=$3
output=$4

## Print to console:
echo "$reference"
echo "$gff"

## Using the input name, create a new variable:
new_name="$(echo "$input" | cut -d '/' -f 2 | cut -d '.' -f 1)"
echo "$new_name"

## Create a tmp directory:
mkdir tmp

## 1) Extract the gene/window coordinates:
echo "Step 1: Extracting the gene/window coordinates"
cut -f 1-3 $input | tail -n +2 | sed 's/\t/-/g' > tmp/"$new_name".bed
cd tmp
echo "Step 1: Extracting the gene/window coordinates - complete"

## 2) Count the number of ambigous bases in each window:
echo "Step 2: Counting number of ambiguous bases"
while read line;
do
echo "$line";
echo "$line" | sed 's/-/\t/g' > "$line".bed;
fastaFromBed -fi ../"$reference" -bed "$line".bed -fo "$line".fasta;
seqtk comp "$line".fasta >> ../results/bed_files_100kb_comp.txt;
done < "$new_name".bed

cd ../

echo "Step 2: Counting number of ambiguous bases - complete!"

## 3) Remove windows where bases make up at least 10 percent as may be the result of gaps/technical artefacts;
Rscript run_filter_windows.R results/bed_files_100kb_comp.txt results/bed_files_100kb_comp_filtered.txt

## 4) Run zscore analysis to identify windows with sig. higher or lower genetic diversity compared to rest of genome:
Rscript run_zscore_analysis.R results/bed_files_100kb_comp_filtered.txt results/"$output"

sort results/"$output"_all_sig.txt | uniq > results/"$output"_all_sig_unique.txt 
sort results/"$output"_all_sig_high.txt | uniq > results/"$output"_all_sig_high_unique.txt
sort results/"$output"_all_sig_low.txt | uniq > results/"$output"_all_sig_low_unique.txt

## 5) Identify genes that overlap with regions of high and low diversity:
## For this step, we can create two new tmp files - one for high div regions and the second for low diversity regions:
mkdir tmp_high
mkdir tmp_low

## Do the same of high diversity regions:
cd tmp_high
cut -f 1 ../results/"$output"_all_sig_high_unique.txt | sed 's/:/-/g' > "$output"_all_sig_high_unique.bed

while read line;
do
echo "$line";
echo "$line" | sed 's/-/\t/g' > "$line".bed;
intersectBed -a ../"$gff" -b "$line".bed | awk '$3=="gene"' | cut -f 9 | cut -d ';' -f 1 | cut -d ':' -f2 > "$line"_gene_list.txt;
done < "$output"_all_sig_high_unique.bed
cd ../

## Do the same of low diversity regions:
cd tmp_low
cut -f 1 ../results/"$output"_all_sig_low_unique.txt | sed 's/:/-/g' > "$output"_all_sig_low_unique.bed

while read line;
do
echo "$line";
echo "$line" | sed 's/-/\t/g' > "$line".bed;
intersectBed -a ../"$gff" -b "$line".bed | awk '$3=="gene"' | cut -f 9 | cut -d ';' -f 1 | cut -d ':' -f2 > "$line"_gene_list.txt;
done < "$output"_all_sig_low_unique.bed
cd ../

## 6) Check if canonical immune genes are found in regions of high and low diversity:
grep -f data/immune_genes.txt tmp_high/*gene_list.txt > results/immune_genes_in_high_diversity_regions.txt
grep -f data/immune_genes.txt tmp_low/*gene_list.txt > results/immune_genes_in_low_diversity_regions.txt
