#!/bin/bash

set -euo pipefail

bam=$1
te=$2
out=$3

bedtools intersect -a $1 -b $te | samtools view -h | awk '{if($0~/^@/){print}else{l=length($10); if(l>=24 && l<=32){print}}}' | samtools view -Sb > $out
