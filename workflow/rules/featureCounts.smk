saf_name = os.path.basename(config["annotations"]).replace(".gff", ".saf")

rule gff_extract:
    input:
        config["annotations"]
    output:
        names = os.path.join(featureCount_dir, "annotations", "names.tsv"),
        saf = os.path.join(featureCount_dir, "annotations", saf_name)
    log:
        os.path.join(featureCount_dir, "log", "gff_extract.log")
    script:
        "../scripts/get_gff_info.py"


rule featureCounts:
    input:
        saf = os.path.join(featureCount_dir, "annotations", saf_name),
        bam = expand(os.path.join(map_out_vb_dir, "{sample}.sorted.bam"), sample = basenames)
    output:
        counts = os.path.join(featureCount_dir, "featureCounts.counts.tsv"),
        outbam = expand(os.path.join(featureCount_dir, "{sample}.sorted.bam.featureCounts.bam"), sample = basenames)
    log:
        os.path.join(featureCount_dir, "log", "featureCounts.log")
    threads: 6
    params:
        " -F SAF -M -O -R BAM --fraction"
    shell:
        "featureCounts {params} -T {threads} -a {input.saf} -o {output.counts} {input.bam}"


rule format_fC_output:
    input:
        os.path.join(featureCount_dir, "featureCounts.counts.tsv"),
        map_out_vb_dir
    output:
        os.path.join(featureCount_dir, "counts.tsv")
    log:
        os.path.join(featureCount_dir, "log", "format_fC_output.log")
    script:
        "../scripts/format_fC_output.py"

# rule format_file_names:
#     input:
#         lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.bam.featureCounts.bam")
#     output:
#         os.path.join(featureCount_dir, "{sample}.featureCounts.bam")
#     log:
#         os.path.join(featureCount_dir, "log", "{sample}.format_file_names.log")
#     script:
#         "../scripts/rename.py"