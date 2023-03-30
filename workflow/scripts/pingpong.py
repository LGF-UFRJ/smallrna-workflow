#!/usr/bin/python3

import pysam
from collections import defaultdict


bamfile = snakemake.input[0]
outfile = snakemake.output[0]
outpairs = snakemake.output[1]

# bamfile = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/pirnas.bam"
# outfile = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/pirnas.pingpong.tsv"
# outpairs = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/pirnas.pingpong.pairs.tsv"


bam = pysam.AlignmentFile(bamfile, "rb")
seen = set()
contig_mappings_pos = defaultdict(list)
contig_mappings_neg = defaultdict(list)
for mapping in bam:
    if not mapping.is_unmapped and mapping.query_name not in seen:
        for tag in mapping.get_tags():
            if tag[0] == "NH" and tag[1] == 1:
                ref_name = mapping.reference_name
                read_id = mapping.query_name
                start = mapping.reference_start
                end = mapping.reference_end
                strand = ""
                if mapping.is_reverse:
                    strand = "-"
                    contig_mappings_neg[ref_name].append({"read_id": read_id, 
                                                        "start": start, 
                                                        "end": end,
                                                        "strand": strand})
                else:
                    strand = "+"
                    contig_mappings_pos[ref_name].append({"read_id": read_id, 
                                                        "start": start, 
                                                        "end": end,
                                                        "strand": strand})
        seen.add(mapping.query_name)

del seen

overlap_count = defaultdict(int)
with open(outpairs, "w") as op:
    for contig in contig_mappings_pos:
        for positive in contig_mappings_pos[contig]:
            for negative in contig_mappings_neg[contig]:
                if positive["start"] < negative["start"] or positive["start"] > negative["end"]:
                    break
                elif negative["start"] <= positive["start"] <= negative["end"]:
                    overlap = negative["end"] - positive["start"]
                    if overlap > 0:
                        print(contig, 
                              positive["start"], 
                              positive["end"], 
                              positive["read_id"], 
                              ".", 
                              positive["strand"],
                              negative["start"], 
                              negative["end"], 
                              negative["read_id"], 
                              ".", 
                              negative["strand"],
                              sep="\t",
                              file=op)
                        overlap_count[overlap] += 1

            

with open(outfile, "w") as out:
    print("overlap_length", "count", sep="\t", file=out)
    for overlap in sorted(overlap_count.keys()):
        print(overlap, overlap_count[overlap], sep="\t", file=out)
