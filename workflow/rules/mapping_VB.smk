
index = config["bowtie_vb_index"]

rule bowtie_vb_genome:
    input:
        lambda wildcards: os.path.join(trim_outdir, wildcards.sample + "_trimmed.fq.gz")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sam")
    params:
        # f"-v 3 -k 5 -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
        f"-v 3 -a -m 50 -S -p 6 --best --strata --chunkmbs 1000 {index}"
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.bowtie.log")
    shell:
        "bowtie {params} {input} > {output}"


rule samtools_vb_convert:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sorted.bam")
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.samtools.log")
    shell:
        "samtools view -Sb {input} | samtools sort > {output} && "
        "samtools index {output} && "
        "rm {input} "

rule NH_tag_vb:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sorted.nh.bam")
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.NH_tag_vb.log")
    shell:
        "scripts/add_nh.sh {input}"


rule get_uniquely_vb:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sorted.uniquely.bam")
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.get_uniquely.log")
    conda:
        "../envs/pysam0.yaml"
    script:
        "../scripts/separate_uniquely.py"


rule index_uniquely_vb:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.uniquely.bam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sorted.uniquely.bam.bai")
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    shell:
        "samtools index {input}"


rule count_mapped_vb:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.sorted.counts.tsv")
    conda:
        "../envs/pysam0.yaml"
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.count.log")
    script:
        "../scripts/count_mapped.py"


# rule sambamba_markdup:
#     input:
#         lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.bam")
#     output:
#         os.path.join(map_out_vb_dir, "{sample}.markdup.bam")
#     shell:
#         "sambamba markdup {input} {output}"

