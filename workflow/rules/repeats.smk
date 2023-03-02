
rule get_genome_unassigned:
    input:
        lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.bam.featureCounts.bam")
    output:
        os.path.join(pirna_dir, "{sample}.genes.unassigned.bam")
    log:
        os.path.join(pirna_dir, "log", "{sample}.genome_unassigned.log")
    conda:
        "../envs/env_pysam.yaml"
    script:
        "../scripts/get_unassigned.py"


rule featureCounts_repeats:
    input:
        gtf = config["repeat_annotation"],
        bam = expand(os.path.join(pirna_dir, "{sample}.genes.unassigned.bam"), sample = basenames)
    output:
        counts = os.path.join(pirna_dir, "repeats.featureCounts.counts.tsv"),
        outbam = expand(os.path.join(pirna_dir, "{sample}.genes.unassigned.bam.featureCounts.bam"), sample = basenames)
    log:
        os.path.join(pirna_dir, "log", "repeats.featureCounts.log")
    threads: 6
    params:
        " -F GTF -M -O -R BAM --fraction -g class_id -t repeat"
    shell:
        "featureCounts {params} -T {threads} -a {input.gtf} -o {output.counts} {input.bam}"

# rule format_file_names_repeats:
#     input:
#         lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.bam.featureCounts.bam")
#     output:
#         os.path.join(featureCount_dir, "{sample}.featureCounts.bam")
#     log:
#         os.path.join(featureCount_dir, "log", "{sample}.format_file_names.log")
#     script:
#         "../scripts/rename.py"