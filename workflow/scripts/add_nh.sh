#!/bin/bash

set -euo pipefail

bam=$1
outfile=$(echo $bam | sed s/bam/nh.bam/)

samtools view -h $1 | awk -v OFS="\t" '{if($0!~/^@/){split($NF, xm, ":"); m=xm[3]; o=0; if(m == 0){o=0}else{o=m-1} $NF="NH:i:"o; print}else{print}}'| samtools view -Sb > $outfile

samtools index $outfile