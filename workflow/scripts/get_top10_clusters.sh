#!/bin/bash


cat $1 | cut -f 1,2,3,5 | head -n 11 | awk 'NR>1{OFS="\t"; print $1,$2,$3,$4,255,"."}' > $2
