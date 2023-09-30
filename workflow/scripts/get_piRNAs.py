import pysam
import re


bamfile = snakemake.input[0]
outfile = snakemake.output[0]

repeat_pattern = re.compile(r"::")
nonrepeat_pattern = re.compile(r"Simple_repeat|rRNA|snRNA|Low_complexity")

bam = pysam.AlignmentFile(bamfile, "rb")
outbam = pysam.AlignmentFile(outfile, "w", template=bam)
for mapping in bam:
    if mapping.is_unmapped == False:
        for tag in mapping.get_tags():
            if tag[0] == "XT":
                if repeat_pattern.search(tag[1]) and not nonrepeat_pattern.search(tag[1]):
                    if 24 <= mapping.query_length <= 32:
                        outbam.write(mapping)
bam.close()
outbam.close()
