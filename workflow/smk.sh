#!/bin/bash

module load snakemake deeptools

set -euo pipefail

snakemake -j 6 --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log" --latency-wait 120 --keep-going

