#!/bin/bash

set -euo pipefail

bamfile=$1

samtools view $bamfile | awk '$3!="*"{x[$1]+=1}END{for(i in x){print i"\t"x[i]}}'

