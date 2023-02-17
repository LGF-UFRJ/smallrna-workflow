from csv import reader
import re

gff = snakemake.input[0]
names = snakemake.output[0]
saf = snakemake.output[1]

with open(gff, "r") as gff_fh, open(names, "w") as names_fh, open(saf, "w") as saf_fh:
    print("GeneID", "Chr", "Start", "End", "Strand", sep="\t", file=saf_fh)
    gff_csv = reader(gff_fh, delimiter="\t")
    for line in gff_csv:
        if not line[0].startswith("#"):
            gene_id = re.search(r"^ID=([A-Z]+[0-9]+);", line[8])
            if gene_id:
                name_re = re.search(r"Name=(.+);d", line[8])
                description = re.search(r"description=(.+);e", line[8])
                if name_re:
                    print(gene_id.group(1), name_re.group(1), sep="\t", file=names_fh)
                else:
                    print(gene_id.group(1), 
                            description.group(1).replace(" ", "_"), 
                            sep="\t", 
                            file=names_fh)
                print(gene_id.group(1), 
                        line[0], 
                        line[3], 
                        line[4], 
                        line[6], 
                        sep="\t", 
                        file=saf_fh)