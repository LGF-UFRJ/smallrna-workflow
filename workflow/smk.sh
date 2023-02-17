#!/bin/bash

module load snakemake deeptools

set -euo pipefail

snakemake -j 6 --use-conda --latency-wait 60 --keep-going --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log"

