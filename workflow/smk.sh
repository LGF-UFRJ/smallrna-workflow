#!/bin/bash

set -euo pipefail

read parameter <<< $@

snakemake --cores 12 $parameter --reason --use-conda --latency-wait 60 --keep-going 

