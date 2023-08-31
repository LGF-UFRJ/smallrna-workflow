#!/bin/bash

set -euo pipefail


echo "pairs="$1
echo "sense="$2
echo "antisense="$3
echo "out_sense="$4
echo "out_antisense="$5

awk '{if(NR==FNR){if($8-$2 == 10){ids[$4]=1; ids[$9]=1}}else{if($0~/^@/){print}else if($1 in ids){print $0}}}' $pairs <(samtools view -h $sense) | samtools view -Sb > $out_sense 

awk '{if(NR==FNR){if($8-$2 == 10){ids[$4]=1; ids[$9]=1}}else{if($0~/^@/){print}else if($1 in ids){print $0}}}' $pairs <(samtools view -h $antisense) | samtools view -Sb > $out_antisense 
