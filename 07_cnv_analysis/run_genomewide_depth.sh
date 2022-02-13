#!/bin/sh

#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH -A mulif005c
#SBATCH -p ProdQ

for name in results/*sorted
do
echo "$name";
new_name="$(echo "$name" | cut -d '.' -f 1 )"
samtools depth "$name" > "$new_name".depth.txt
done
