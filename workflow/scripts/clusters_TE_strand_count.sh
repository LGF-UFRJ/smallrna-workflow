#!/bin/bash

set -euo pipefail

te_annotation=$1
clusters=$2

bedtools intersect -wo -a $te_annotation -b $clusters | awk -v OFS="\t" 'BEGIN{print "cluster", "pos", "neg"}{cluster=$10; strand=$6; all_cls[$10]=1; if(strand == "+"){psc[cluster]+=1}else{nsc[cluster]+=1}}END{for(i in all_cls){if(psc[i]==""){psc[i]=0}else if(nsc[i]==""){nsc[i]=0};print i, psc[i], nsc[i]*-1}}' 
