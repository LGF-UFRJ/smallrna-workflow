import subprocess as sp


bamfile = snakemake.input[0]
genome_windows = snakemake.input[1]
clusters = snakemake.output[0]

# bamtobed = f"bamToBed -i {bamfile} | awk '$3-$2>=24 && $3-$2<=30' > {bed}"
# sp.run(bamtobed, shell=True)

awk = "awk '$3-$2>5000{print $0\"\\t\"$3-$2}'"

# get_clusters = f"bedtools intersect -c -a {genome_windows} -b {bed} | awk '$4>=5' | bedtools merge -d 20000 | bedtools intersect -c -a - -b {bed} | {awk} > {clusters}"

# get_clusters = f"bedtools intersect -c -a {genome_windows} -b {bed} | awk '$4>5' | bedtools merge | awk '$3-$2>=5000' | bedtools merge -d 20000 | {awk} > {clusters}"

# get_clusters = f"bedtools intersect -c -a {genome_windows} -b {bed} | awk '$4>25' | bedtools merge -d 20000 | {awk} > {clusters}"

get_clusters = f"samtools view -b {bamfile} | bedtools intersect -c -a {genome_windows} -b - | awk '$4>25' | bedtools merge -d 1 | {awk} > {clusters}"

sp.run(get_clusters, shell=True)
