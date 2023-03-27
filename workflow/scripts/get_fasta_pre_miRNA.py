import csv
import os
import subprocess as sp

annot = snakemake.input[0]
genome = snakemake.input[1]

output = snakemake.output[0]

tmp = "./temp.gtf"

with open(annot, "r") as annot_fh, open(tmp, "w") as tmp_fh:
    annot_csv = csv.reader(annot_fh, delimiter="\t")
    for line in annot_csv:
        if not line[0].startswith("#") and line[2] == "pre_miRNA":
            print(*line, sep="\t", file=tmp_fh)

module = "module load bedtools2 && "
command = f"bedtools getfasta -s -fi {genome} -fo {output} -bed {tmp}"
sp.run(module + command, shell=True)
os.remove(tmp)