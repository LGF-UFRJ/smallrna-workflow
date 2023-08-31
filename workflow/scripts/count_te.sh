#!/bin/bash

set -euo pipefail

mapreadcount=$1
featurecount=$2

awk -v OFS="\t" 'BEGIN{print "repeat", "total", "pos", "neg"}{if(NR==FNR){x[$1]=$2}else{if(length($10)>=24 && length($10)<=32){split($NF, r, "::"); if(x[$1] != ""){c=(1/x[$1]); te_count[r[2]]+=c; if($2==16){te_count_n[r[2]]+=c}else{te_count_p[r[2]]+=c}}}}}END{for(i in te_count){if(i != ""){print i, te_count[i], te_count_p[i], te_count_n[i]}}}' <(cat $mapreadcount) <(samtools view $featurecount)


