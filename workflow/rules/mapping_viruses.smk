
index_vir = config["bowtie_virus_index"]

rule get_unmapped:
    input:
        os.path.join(map_out_dir, "{sample}.sorted.bam")
    output:
        os.path.join(map_out_virus_dir, "{sample}.sorted.unmapped.fq")
    shell:
        "samtools view -f 4 {input} | samtools fastq > {output}"


rule bowtie_viruses:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".sorted.unmapped.fq")
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.sam")
    threads: 6
    params:
        # f"-v 3 -k 5 -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
        f"-v 3 -a -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index_vir}"

    log:
        os.path.join(map_out_virus_dir, "log", "{sample}.bowtie.log")
    shell:
        "bowtie {params} {input} > {output}"


rule samtools_viruses_convert:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sam")
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.bam")
    log:
        os.path.join(map_out_virus_dir, "log", "{sample}.samtools.log")
    shell:
        "samtools view -Sb {input} | samtools sort > {output} && "
        "samtools index {output} && "
        "rm {input} "


# rule index_uniquely_hic:
#     input:
#         lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam")
#     output:
#         os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam.bai")
#     log:
#         os.path.join(map_out_dir, "log", "{sample}.index_uniquely.log")
#     shell:
#         "samtools index {input}"


# rule count_mapped_viruses:
#     input:
#         lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam")
#     output:
#         os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.counts.tsv")
#     conda:
#         "../envs/pysam0.yaml"
#     log:
#         os.path.join(map_out_virus_dir, "log", "{sample}.count.log")
#     script:
#         "../scripts/count_mapped.py"

rule length_distribution_viruses:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
        os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.counts.tsv")
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.lengths.tsv")
    params:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
    conda:
        "../envs/pysam0.yaml"
    log:
        os.path.join(map_out_virus_dir, "log", "{sample}.lengths.log")
    script:
        "../scripts/length_distribution.py"

rule plot_length_distribution_viruses:
    input:
        expand(os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.lengths.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(map_out_virus_dir, "ridges.lengths.plot.png"),
        os.path.join(map_out_virus_dir, "line.lengths.plot.png"),
    log:
        os.path.join(map_out_virus_dir, "log", "lengths.plot.log")
    script:
        "../scripts/plot_length_distribution.R"


# # rule sambamba_markdup:
# #     input:
# #         lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam")
# #     output:
# #         os.path.join(map_out_dir, "{sample}.markdup.bam")
# #     shell:
# #         "sambamba markdup {input} {output}"

