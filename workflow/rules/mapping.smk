
index = config["bowtie_index"]

rule bowtie_genome:
    input:
        lambda wildcards: os.path.join(trim_outdir, wildcards.sample + "_trimmed.fq.gz")
    output:
        os.path.join(map_out_dir, "{sample}.sam")
    params:
        f"-v 3 -k 5 -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
    log:
        os.path.join(map_out_dir, "log", "{sample}.bowtie.log")
    shell:
        "bowtie {params} {input} > {output}"


rule samtools_convert:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.bam")
    log:
        os.path.join(map_out_dir, "log", "{sample}.samtools.log")
    shell:
        "samtools view -Sb {input} | samtools sort > {output} && "
        "samtools index {output} && "
        "rm {input} "



rule count_mapped:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
    conda:
        "../envs/small.yaml"
    log:
        os.path.join(map_out_dir, "log", "{sample}.count.log")
    script:
        "../scripts/count_mapped.py"


# rule sambamba_markdup:
#     input:
#         lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
#     output:
#         os.path.join(map_out_dir, "{sample}.markdup.bam")
#     shell:
#         "sambamba markdup {input} {output}"

