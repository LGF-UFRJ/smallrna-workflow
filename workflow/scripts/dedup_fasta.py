import sys

seq_dict = {}

read_id = ""
seq = ""
for line in sys.stdin:
    if line.startswith(">"):
        if read_id == "":
            read_id = line.strip()
        else:
            seq_dict[read_id] = seq
            read_id = line.strip()
            seq = ""
    else:
        seq += line.strip()

seq_dict[read_id] = seq

for read_id in seq_dict:
    print(read_id, seq_dict[read_id], sep="\n")