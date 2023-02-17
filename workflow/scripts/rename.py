import os

if not os.path.exists(snakemake.output[0]):
    os.rename(snakemake.input[0], snakemake.output[0])