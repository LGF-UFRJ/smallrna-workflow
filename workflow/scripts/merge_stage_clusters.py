import subprocess as sp
import re
import os.path

clusters = snakemake.input

output = snakemake.output

for out in output:
    stage = os.path.basename(out.replace(".clusters.nomir.bed", ""))
    pat1 = re.compile(fr"{stage}1")
    pat2 = re.compile(fr"{stage}2")
    cl1 = [i for i in clusters if pat1.search(i)][0]
    cl2 = [i for i in clusters if pat2.search(i)][0]
    command = f"cat {cl1} {cl2} | bedtools sort | bedtools merge > {out}"
    sp.run(command, shell=True)
