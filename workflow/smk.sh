#!/bin/bash

module load slurm snakemake/5.4.3  deeptools R trim_galore bowtie samtools MultiQC subread bedtools2 UCSCtools

set -euo pipefail

snakemake -j 8 --reason --use-conda --latency-wait 60 --keep-going --cluster "SlurmEasy -n {rule} -t {threads} -l cluster_log -x deep17"

