map_out_dir = os.path.join(config["results"], "mapping/")
index = config["bowtie_index"]



rule bowtie_genome:
    input:
        lambda wildcards: os.path.join(trim_outdir, wildcards.sample + "_trimmed.fq.gz")
    output:
        os.path.join(map_out_dir, "{sample}.sam")
    params:
        f"-v 3 -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
    log:
        os.path.join(map_out_dir, "log", "{sample}.bowtie.log.txt")
    shell:
        "bowtie {params} {input} 1> {output} 2> {log}"

rule samtools_convert:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.bam")
    shell:
        "samtools view -Sb {input} | samtools sort > {output} && "
        "samtools index {output} && "
        "rm {input} "
