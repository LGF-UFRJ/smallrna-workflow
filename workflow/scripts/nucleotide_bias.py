import pysam

bamfile = snakemake.input[0]
outfile = snakemake.output[0]

A, C, T, G, N = {}, {}, {}, {}, {} 
total = 0
for i in range(1, 33):
    A[i] = 0
    C[i] = 0
    T[i] = 0
    G[i] = 0
    N[i] = 0
bam = pysam.AlignmentFile(bamfile, "rb")
seen = set()
for mapping in bam:
    if not mapping.is_unmapped and mapping.query_name not in seen:
        for position, nucleotide in enumerate(mapping.get_forward_sequence()):
            position += 1
            if nucleotide == "A":
                A[position] += 1
            elif nucleotide == "T":
                T[position] += 1
            elif nucleotide == "G":
                G[position] += 1
            elif nucleotide == "C":
                C[position] += 1
            elif nucleotide == "N":
                N[position] += 1
        total += 1
        seen.add(mapping.query_name)


def pct(x, total): return (x*100)/total
with open(outfile, "w") as out:
    print("position", "A", "C", "T", "G", "N", sep="\t", file=out)
    for position in A:
        print(position,
                pct(A[position], total),
                pct(C[position], total),
                pct(T[position], total), 
                pct(G[position], total),
                pct(N[position], total),
                sep="\t",
                file=out)