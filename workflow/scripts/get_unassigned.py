
# I'm reading the file 2 times here because I want to get reads that are 
# "truly unassigned". So, multimapped reads with at least one mapping in an 
# annotated feature are consider assigned and discarded. For that, I need first
# to map the ones considered "assigned" and then use that to filter the 
# "unassigned" ones. Maybe there is a better and less complex way to do this, 
# but I cannot think of it right now  ¯\_(ツ)_/¯

import pysam

assigned = set()
bam = pysam.AlignmentFile(snakemake.input[0], "rb")
outbam = pysam.AlignmentFile(snakemake.output[0], "w", template=bam)

for mapping in bam:
    if mapping.is_unmapped == False:
        for tag in mapping.get_tags():
            if tag == "XS:Z:Assigned":
                assigned.add(mapping.query_name)

for mapping in bam:
    if mapping.is_unmapped == False and mapping.query_name not in assigned:
        outbam.write(mapping)

