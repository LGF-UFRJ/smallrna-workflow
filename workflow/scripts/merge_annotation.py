from csv import reader
import subprocess as sp

repeats_gtf = snakemake.input[0]
features_saf = snakemake.input[1]
repeats_saf = snakemake.output[0]
merged_saf = snakemake.output[1]

with open(repeats_gtf, "r") as gtf, open(repeats_saf, "w") as saf:
    print("GeneID", "Chr", "Start", "End", "Strand", sep="\t", file=saf)
    gtf_read = reader(gtf, delimiter="\t")
    for line in gtf_read:
        info = line[8].split()
        geneid = info[1].replace("\"", "").replace(";", "")
        print(geneid, line[0], line[3], line[4], line[6], sep="\t", file=saf)

sp.run(f"cat {features_saf} {repeats_saf} > {merged_saf}", shell=True)