import pandas as pd
import subprocess as sp
import os.path

pd.set_option('display.max_colwidth', None)
samplesheet = pd.read_csv(snakemake.input[0], sep="\t")
outdir = snakemake.params[0]
sp.run(f"mkdir -p {outdir}", shell=True)

for name in samplesheet["name"]:
    path = samplesheet[samplesheet["name"] == name]["path"].to_string(index = False)
    outname = os.path.join(outdir, name + ".fastq.gz")
    command = f"ln -s {path} {outname}"
    run = sp.run(command, shell=True, capture_output=True)
    print(command)
    print(run.stdout)
    print(run.stderr)
