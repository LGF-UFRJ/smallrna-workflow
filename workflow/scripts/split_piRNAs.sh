#!/bin/bash

set -euo pipefail


pairs=$1
sense=$2
antisense=$3
out_sense=$4
out_antisense=$5

awk '{if(NR==FNR){ids[$4]=1; ids[$9]=1}else{if($0~/^@/){print}else if($1 in ids){print $0}}}' $pairs <(samtools view -h $sense) | samtools view -Sb > $out_sense 

awk '{if(NR==FNR){ids[$4]=1; ids[$9]=1}else{if($0~/^@/){print}else if($1 in ids){print $0}}}' $pairs <(samtools view -h $antisense) | samtools view -Sb > $out_antisense 
