#!/bin/bash
set -euo pipefail

samtools view $1 | awk 'length($10)>=24 && length($10)<=32{x[$1]+=1}END{for(i in x){if(x[i] == 1){print i}}}' > $2
