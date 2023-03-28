#!/usr/bin/python3

from csv import reader
import os.path

fc_out = snakemake.input[0] 
# map_path = snakemake.params[0]
fmt_out = snakemake.output[0]   

with open(fc_out, "r") as fc, open(fmt_out, "w") as out:
    fc_csv = reader(fc, delimiter="\t")
    for line in fc_csv:
        if not line[0].startswith("#"):
            if line[0].startswith("Geneid"):
                fmt_header = [os.path.basename(i).replace(".sorted.bam", "") for i in line]
                print(*fmt_header, sep="\t", file=out)
            else:
                print(*line, sep="\t", file=out)

