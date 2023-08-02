import subprocess as sp
import pandas as pd
import argparse
import re

# Normalization function
def normalize(n, total):
    return (n * 1_000_000) / total

# Arguments
parser = argparse.ArgumentParser(description="")
parser.add_argument("-w", "--window-file", help="File with counts for each window")
parser.add_argument("-o", "--output", help="Output file with clusters")
parser.add_argument("-c", "--count-files", nargs = "*", help="files with counts for each sample")
args = vars(parser.parse_args())

# Inputs
window_file = args["window_file"]
output_file = args["output"]
count_files = args["count_files"]


# Parsing total mapped read counts for each sample
total = {}
sample_pat = re.compile(r"([a-z]+[0-9]).sorted")
for f in count_files:
    sample = sample_pat.search(f).group(1)
    countsheet = pd.read_csv(f, sep="\t", names = ["category", "count"])
    total[sample] = int(countsheet[countsheet["category"] == "total"]["count"])

# Formatting header of the window count file
windowsheet = pd.read_csv(window_file, sep="\t")
columns = windowsheet.columns
rename_col = {} 
for i in columns:
    noq = re.sub(r"\'", "", i)
    noh = re.sub(r"#", "", noq)
    noe = re.sub(r"\.piRNAs\.uniquely\.bam", "", noh)
    rename_col[i] = noe
windowsheet_fmt = windowsheet.rename(columns=rename_col)
columns_fmt = windowsheet_fmt.columns
    
# Normalizing window count values
for i in range(3, len(columns_fmt)):
    sample = columns_fmt[i]
    sample_vals = windowsheet_fmt[sample]
    sample_norm = sample_vals.apply(normalize, args = (total[sample], ))
    windowsheet_fmt[sample] = sample_norm

# Saving normalized counts
windowsheet_fmt.to_csv(output_file, sep="\t", index=False)
