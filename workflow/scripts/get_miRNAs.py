import pysam
import pandas as pd
import re

bamfile = snakemake.input[0]
names = snakemake.input[1]
outfile = snakemake.output[0]


names_df = pd.read_csv(names, sep="\t", header=None)
mir_pattern = re.compile(r"mirna|mir", re.I)
mir_filter = [True if mir_pattern.search(i) else False for i in names_df[1]]
mir_ids = names_df[mir_filter][0].to_list()

bam = pysam.AlignmentFile(bamfile, "rb")
outbam = pysam.AlignmentFile(outfile, "wb", template=bam)
for mapping in bam:
    if mapping.is_unmapped == False:
        for tag in mapping.get_tags():
            if tag[0] == "XT" and tag[1] in mir_ids:
                outbam.write(mapping)
bam.close()
