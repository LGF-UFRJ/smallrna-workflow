#!/bin/bash

set -euo pipefail

read parameter <<< $@

nice -n 10 snakemake --cores 24 $parameter --reason --use-conda --latency-wait 120 --keep-going 

