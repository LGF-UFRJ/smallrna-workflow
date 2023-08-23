#!/bin/bash

set -euo pipefail

top10=$1
outdir=$2

for i in $(cut -f 4 $top10 | sort | uniq); do cat $top10 | awk -v i=$i '$4==i{split($10, x, "::"); te_count[x[2]] += 1}END{for(i in te_count){OFS="\t"; print i, te_count[i]}}' > "${outdir}/${i}.TEs.count.tsv"; done 
