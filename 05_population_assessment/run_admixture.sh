#!/bin/sh

#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

for K in {1..20};
do
admixture --cv results/out_irish_hetremoved_raresnpsremoved.ann_plink.bed \
-j20 $K | tee results/log${K}.out;
done
