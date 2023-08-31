#!/bin/bash

set -euo pipefail

bamfile=$1

samtools view -h $bamfile | awk '{if($0~/^@/){print}else{if(x[$1]!=1){print; x[$1]=1}}}' | samtools fasta -   


