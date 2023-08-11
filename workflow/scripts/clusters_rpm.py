import pandas as pd
import glob
import os.path

def get_total(tsv, total_dict):
    basename = os.path.basename(tsv).replace(".sorted.counts.tsv", "")
    with open(tsv, "r") as t:
        for line in t:
            col = line.strip().split("\t")
            if col[0].startswith("total"):
                total_dict[basename] = int(col[1])
    return total_dict

def normalize(x, total, factor=1_000_000):
    return (x*factor) / total

table_file = snakemake.input[0]
table_pos_file = snakemake.input[1]
table_neg_file = snakemake.input[2]
counts = snakemake.input[3:]
outtable = snakemake.output[0]

table = pd.read_csv(table_file, sep="\t")
table_pos = pd.read_csv(table_pos_file, sep="\t")
table_neg = pd.read_csv(table_neg_file, sep="\t")

header = [i.replace("'", "").replace("#", "") for i in table.columns]
table.columns = header
table_pos.columns = header
table_neg.columns = header

total_dict = {}
for c in counts:
    total_dict = get_total(c, total_dict)

bed = ["chr", "start", "end"]
for colname in table.columns:
    if colname not in bed:
        sample = colname.replace(".sorted.bam", "")
        total = total_dict[sample]
        norm_list = [normalize(i, total) for i in table[colname]]
        table[colname] = norm_list
        # norm_pos_list = [normalize(i, total, factor=100) for i in table_pos[colname]]
        # table_pos[colname] = norm_pos_list

lengths = table["end"] - table["start"]
table.insert(3, "length", lengths)

sample_matrix = table.iloc[:,4:]
sample_mean = sample_matrix.mean(axis=1)
table.insert(4, "mean_rpm", sample_mean)

sample_pos_matrix = table_pos.iloc[:,3:]
sample_neg_matrix = table_neg.iloc[:,3:]
sample_pos_mean = sample_pos_matrix.mean(axis=1)
sample_neg_mean = sample_neg_matrix.mean(axis=1)
table_pos["mean"] = sample_pos_mean
table_neg["mean"] = sample_neg_mean

table_pos_neg = table_pos.merge(table_neg, how="left", on = ["chr", "start", "end"])
total = table_pos_neg["mean_x"] + table_pos_neg["mean_y"]
pos_ratio = table_pos_neg["mean_x"] / total

cl_type = []
for i in pos_ratio:
    if i >= 0.8:
        cl_type.append("uni:pos")
    elif i <= 0.2:
        cl_type.append("uni:neg")
    else:
        cl_type.append("dual")

table_pos_neg["type"] = cl_type

final = table.merge(table_pos_neg[["chr", "start", "end", "type"]], 
                    how="left", 
                    on = ["chr", "start", "end"])
# final["factor"] = final["mean_rpm"] * final["length"] 
final = final.sort_values(by = ["length"], ascending=False)
names = ["RPCL"+str(i+1) for i in range(len(final))]
final.insert(5, "name", names)
final = final[["chr", "start", "end", "length", "name", "mean_rpm", "type"]]
final.to_csv(outtable, index=False, sep="\t")
