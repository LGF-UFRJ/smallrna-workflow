import pysam
import pandas as pd

def normalize(x, total, factor=1_000_000):
    return (x/total) * factor

bamfile = snakemake.input[0]
outtable = snakemake.output[0]
counts = snakemake.params[0]

counts_df = pd.read_csv(counts, sep = "\t", names = ["type", "freq"])
nmapped = int(counts_df[counts_df["type"] == "total"]["freq"])

bam = pysam.AlignmentFile(bamfile, "rb")
seen = set()
zeros = [0] * len([*range(18, 76)])
len_dist = {"length": [*range(18, 76)],
            "total": zeros.copy(),
            "pos": zeros.copy(),
            "neg": zeros.copy()}
for mapping in bam:
    if not mapping.is_unmapped and mapping.query_name not in seen:
        position = mapping.query_length - 18
        len_dist["total"][position] += 1
        if not mapping.is_reverse:
            len_dist["pos"][position] += 1
        else:
            len_dist["neg"][position] -= 1
        seen.add(mapping.query_name)

len_dist["total"] = [normalize(i, nmapped) for i in len_dist["total"]]
len_dist["pos"] = [normalize(i, nmapped) for i in len_dist["pos"]]
len_dist["neg"] = [normalize(i, nmapped) for i in len_dist["neg"]]
outdf = pd.DataFrame(len_dist)

outdf.to_csv(outtable, sep = "\t", index = False)