#!/bin/bash

set -euo pipefail

cat $1 | awk -v OFS="\t" '{print $1,$4-1,$5,substr($10,2,length($10)-3), 255, $7}' 
