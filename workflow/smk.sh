#!/bin/bash

module load slurm snakemake/7.0.4 deeptools R

set -euo pipefail

snakemake -j 6 --use-conda --latency-wait 60 --keep-going --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log"

