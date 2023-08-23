#!/bin/bash

set -euo pipefail

bam=$1
window=$2
basename=$(echo $bam | sed -e 's/\.piRNAs.uniquely.bam//' -e 's/clusters\/piRNAs\//clusters\//')

merged=${basename}.merged.bed
wcount=${basename}.window_count.tsv
clusters=${basename}.clusters.bed

pirnaDensity=0.12
minClusterSize=5000
maxClusterDist=20000

echo ""
echo -e "\tMerging piRNA reads..."
# bedtools merge -s -c 1 -o collapse -i $bam | awk '{split($4, x, ","); if(length(x) > 1){print}}' | grep -v "\-1" > $merged
bedtools merge -s -c 1 -o collapse -i $bam | grep -v "\-1" > $merged

echo -e "\tCalculating coverage of sliding window..."
bedtools coverage -a $window -b $merged > $wcount

echo -e "\tGetting clusters..."
cat $wcount |\
    awk -v pdens=$pirnaDensity '$7 > pdens' |\
    bedtools merge |\
    awk -v mincl=$minClusterSize '$3-$2 > mincl' |\
    bedtools merge -d $maxClusterDist |\
    awk '{print $0"\t"$3-$2}' > \
    $clusters

echo -e "\tDONE!"
echo ""
echo -e "\tResults in file: $clusters"
echo ""