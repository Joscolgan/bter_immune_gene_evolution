
## A script to take a tab-delimited text file as input:
## - Extract the gene/window coordinates;
## - Count the number of ambigious bases in each window;
## - Remove windows where bases make up at least 10 percent as may be the result of gaps/technical artefacts;
## - Run zscore analysis to identify windows with sig. higher or lower genetic diversity compared to rest of genome:

## Take input and output as arguments from the command line:
input=$1
reference=$2
output=$3

## 1) Extract the gene/window coordinates:
cut -f 1-3 $input | sed 's/\t/-/g' > "$input".bed

## 2) Count the number of ambigious bases in each window:
while read line;
do
echo "$line" | sed 's/-/\t/g' > "$line".bed
fastaFromBed -fi $reference -bed "$line".bed -fo "$name".fasta;
seqtk comp "$name".fasta >> results/bed_files_100kb_comp.txt ; done
done < "$input".bed

## 3) Remove windows where bases make up at least 10 percent as may be the result of gaps/technical artefacts;
Rscript filter_windows.R results/bed_files_100kb_comp.txt results/bed_files_100kb_comp_filtered.txt

## 4) Run zscore analysis to identify windows with sig. higher or lower genetic diversity compared to rest of genome:
Rscript run_zscore_analysis.R results/bed_files_100kb_comp_filtered.txt results/"$output"

