import subprocess as sp
import pandas as pd
import argparse
import re

count_files = [
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/pvs1.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/pvs2.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/egg1.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/egg2.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/emb1.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/emb2.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/nym1.sorted.counts.tsv",
    "/home5/attilio/Tarcisio/smallrna_rhodnius_2023/results/mapping/nym2.sorted.counts.tsv",
]

# Normalization function
def normalize(n, total):
    return (n * 1_000_000) / total

# Arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("-o", "--output", help="Output file")
parser.add_argument("-c", "--count-files", nargs = "*", help="files with counts for each sample")
args = vars(parser.parse_args())

# Inputs
output_file = args["output"]
count_files = args["count_files"]

# Parsing total mapped read counts for each sample
scallling_factors = {}
sample_pat = re.compile(r"([a-z]+[0-9]).sorted")
for f in count_files:
    sample = sample_pat.search(f).group(1)
    countsheet = pd.read_csv(f, sep="\t", names = ["category", "count"])
    total = int(countsheet[countsheet["category"] == "total"]["count"])
    scallling_factors[sample] = normalize(1, total)

output = pd.DataFrame({"sample": scallling_factors.keys(), "scallingFcator": scallling_factors.values()})

output.to_csv(output_file, index=False, sep="\t")
