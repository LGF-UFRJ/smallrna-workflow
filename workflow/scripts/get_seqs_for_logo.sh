#!/bin/bash

set -euo pipefail

fasta=$1

cat $fasta |\
  grep -v ">" |\
  awk '{seq=substr($0, 1, 20); split(seq,x,""); for(i=1; i<=length(x); i++){if(x[i] == "A" || x[i] == "T"){ta+=1}}; if(ta/length(x) < 0.8){print seq}; ta = 0}' |\
  tr 'T' 'U' 

#awk '{print substr($0, 1, 20)}' |\
