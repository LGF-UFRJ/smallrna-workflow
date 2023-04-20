#!/bin/bash

set -euo pipefail

final=$1
index=$2
sorted=$3
output=$4


cat $final | awk 'NR>1{OFS="\t"; print $1,$2,$3,$5,255,"."}' |  bedtools sort > $sorted
bedToBigBed $sorted $index $output