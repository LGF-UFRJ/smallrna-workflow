#!/bin/bash

set -euo pipefail

input=$1

cl_files=$(realpath "${input}/"* | grep "[12].clusters.nomir.bed")

for i in $cl_files; 
do
    for j in $(echo $cl_files | tr ' ' '\n' | grep -v $i); 
    do 
        bedtools intersect -a $i -b $j | wc -l | awk -v i=${i%%.*} -v j=${j%%.*} '{split(i, fi, "/"); split(j, fj, "/"); li=length(fi); lj=length(fj); print fi[li]"_"fj[lj]"\t"$1}'
    done
done