#!/bin/bash

set -euo pipefail

NF=$(head -n 1 $1 | awk '{print NF-1}')
header=($(head -n 1 $1 | sed "s/'//g"))
outdir=$2
threshold=1
distance=5000

for i in $(seq 3 $NF)
do
    field=$(bc -l <<< "$i + 1")
    sample=${header[i]%%.*}
    sampleout="${outdir}/${sample}.clusters.bed"
    cat $1 | awk -v f=$field -v t=$threshold 'NR>1 && $(f)>t{OFS="\t"; print $1,$2,$3}' | bedtools sort | bedtools merge -d 1 | awk '$3-$2>5000' | bedtools merge -d $distance | awk '$3-$2>5000{print $0"\t"$3-$2}' > $sampleout
done
