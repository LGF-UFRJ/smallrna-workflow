
# I'm reading the file 2 times here because I want to get reads that are 
# "truly unassigned". So, multimapped reads with at least one mapping in an 
# annotated feature are consider assigned and discarded. For that, I need first
# to map the ones considered "assigned" and then use that to filter the 
# "unassigned" ones. Maybe there is a better and less complex way to do this, 
# but I cannot think of it right now  ¯\_(ツ)_/¯

import subprocess; subprocess.run("echo $PATH", shell=True)
import pysam

assigned = set()

bam = pysam.AlignmentFile(snakemake.input[0], "rb")
for mapping in bam:
    if mapping.is_unmapped == False:
        for tag in mapping.get_tags():
            if tag == "XS:Z:Assigned":
                assigned.add(mapping.query_name)
bam.close()

bam = pysam.AlignmentFile(snakemake.input[0], "rb")
outbam = pysam.AlignmentFile(snakemake.output[0], "wb", template=bam)
for mapping in bam:
    if mapping.is_unmapped == False and mapping.query_name not in assigned:
        outbam.write(mapping)

bam.close()
outbam.close()