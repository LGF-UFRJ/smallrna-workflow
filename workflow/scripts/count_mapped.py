import pysam


counted = set()
total_mapped = 0
unique = 0
multi = 0
bam = pysam.AlignmentFile(snakemake.input[0], "rb")

for mapping in bam:
    if mapping.is_unmapped == False and mapping.query_name not in counted:
        total_mapped += 1
        # xm = [i for i in mapping.get_tags() if i[0] == "XM"][0]
        xm = mapping.get_tags()[-1]
        # breakpoint()
        if xm[1] == 2:
            unique += 1
        else:
            multi += 1
            counted.add(mapping.query_name)

with open(snakemake.output[0], "w") as outfile:
    print(f"total\t{total_mapped}", file=outfile)
    print(f"unique\t{unique}", file=outfile)
    print(f"multi\t{multi}", file=outfile)
