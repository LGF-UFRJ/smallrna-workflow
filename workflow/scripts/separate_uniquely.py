import pysam

baminput = snakemake.input[0]
bamoutput = snakemake.output[0]

bam = pysam.AlignmentFile(baminput, "rb")
bamunique = pysam.AlignmentFile(bamoutput, "wb", template=bam)
for mapping in bam:
    if mapping.is_unmapped == False:
        for tag in mapping.get_tags():
            if tag[0] == "XM" and tag[1] == 2:
                bamunique.write(mapping)
bam.close()
bamunique.close()