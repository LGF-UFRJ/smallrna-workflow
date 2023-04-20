#!/bin/bash

set -euo pipefail

samtools view -h $1 | awk '{if($1!~/^@/){if(length($10)>=24 && length($10)<=32){print}}else{print}}' | samtools view -Sb > $2
samtools index $2