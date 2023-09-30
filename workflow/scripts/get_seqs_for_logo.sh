#!/bin/bash

set -euo pipefail

fasta=$1

if [ -z $(cat $fasta) ]; then
  echo ""
else
  cat $fasta |\
    grep -v ">" |\
    awk '{if(length($0) < 20){next}; seq=substr($0, 1, 20); split(seq,x,""); for(i=1; i<=length(x); i++){if(x[i] == "A" || x[i] == "T"){ta+=1}}; if(ta/length(x) < 0.8){print seq}; ta = 0}' |\
    tr 'T' 'U' 
fi

#awk '{print substr($0, 1, 20)}' |\
