genome = snakemake.input[0]
out = snakemake.output[0]

with open(genome, "r") as fasta, open(out, "w") as out_fh:
    for line in fasta:
        if line.startswith(">"):
            header = line.split()
            print(header[0], file=out_fh)
        else:
            print(line.strip(), file=out_fh)