
index = config["bowtie_hic_index"]

rule bowtie_hic_genome:
    input:
        lambda wildcards: os.path.join(trim_outdir, wildcards.sample + "_trimmed.fq.gz")
    output:
        os.path.join(map_out_dir, "{sample}.sam")
    threads: 6
    params:
        # f"-v 3 -k 5 -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
        f"-v 3 -a -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"

    log:
        os.path.join(map_out_dir, "log", "{sample}.bowtie.log")
    shell:
        "bowtie {params} {input} > {output}"


rule samtools_hic_convert:
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


rule get_uniquely_hic:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam")
    log:
        os.path.join(map_out_dir, "log", "{sample}.get_uniquely.log")
    script:
        "../scripts/separate_uniquely.py"


rule index_uniquely_hic:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam.bai")
    log:
        os.path.join(map_out_dir, "log", "{sample}.index_uniquely.log")
    shell:
        "samtools index {input}"


rule count_mapped_hic:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
    log:
        os.path.join(map_out_dir, "log", "{sample}.count.log")
    script:
        "../scripts/count_mapped.py"

rule length_distribution_hic:
    input:
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_dir, "{sample}.sorted.lengths.tsv")
    params:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
    log:
        os.path.join(map_out_dir, "log", "{sample}.lengths.log")
    script:
        "../scripts/length_distribution.py"

rule plot_length_distribution_hic:
    input:
        expand(os.path.join(map_out_dir, "{sample}.sorted.lengths.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(map_out_dir, "ridges.lengths.plot.png"),
        os.path.join(map_out_dir, "line.lengths.plot.png"),
    log:
        os.path.join(map_out_dir, "log", "lengths.plot.log")
    script:
        "../scripts/plot_length_distribution.R"


# rule sambamba_markdup:
#     input:
#         lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
#     output:
#         os.path.join(map_out_dir, "{sample}.markdup.bam")
#     shell:
#         "sambamba markdup {input} {output}"

