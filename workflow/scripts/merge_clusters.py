import subprocess as sp

concat = snakemake.input
clusters = snakemake.output[0]

awk = "awk '{OFS=\"\\t\"; print $0,$3-$2,\".\",\".\"}'"
sp.run(f"cat {concat} | bedtools sort | bedtools merge | {awk} | sort -rnk4 > {clusters}", shell=True)

            
