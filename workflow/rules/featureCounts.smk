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


rule merge_annotation:
    input: 
        repeat_gtf = config["repeat_annotation"],
        features_saf = os.path.join(featureCount_dir, "annotations", saf_name),
    output: 
        repeat_saf = os.path.join(featureCount_dir, "annotations", "repeats.saf"),
        merged_saf = os.path.join(featureCount_dir, "annotations", "merged_annotation.saf")
    script:
        "../scripts/merge_annotation.py" 


rule featureCounts:
    input:
        # saf = os.path.join(featureCount_dir, "annotations", saf_name),
        saf = os.path.join(featureCount_dir, "annotations", "merged_annotation.saf"),
        # bam = expand(os.path.join(map_out_vb_dir, "{sample}.sorted.bam"), sample = samplesheet["name"])
        bam = expand(os.path.join(map_out_vb_dir, "{sample}.sorted.nh.bam"), sample = samplesheet["name"])
    output:
        counts = os.path.join(featureCount_dir, "featureCounts.counts.tsv"),
        outbam = expand(os.path.join(featureCount_dir, "{sample}.sorted.nh.bam.featureCounts.bam"), sample = samplesheet["name"])
    log:
        os.path.join(featureCount_dir, "log", "featureCounts.log")
    threads: 6
    params:
        " -F SAF -M -O -R BAM --fraction"
    shell:
        "featureCounts {params} -T {threads} -a {input.saf} -o {output.counts} {input.bam}"

rule format_fC_output:
    input:
        os.path.join(featureCount_dir, "featureCounts.counts.tsv")
    output:
        os.path.join(featureCount_dir, "counts.tsv")
    log:
        os.path.join(featureCount_dir, "log", "format_fC_output.log")
    script:
        "../scripts/format_fC_output.py"


rule plot_library_profile:
    input:
        expand(os.path.join(map_out_vb_dir, "{sample}.sorted.counts.tsv"), sample = samplesheet["name"]),
    params:
        config["annotations"],
        os.path.join(featureCount_dir, "counts.tsv"),
    output:
        os.path.join(featureCount_dir, "libraries_profile.plot.png")
    log:
        os.path.join(featureCount_dir, "log", "libraries_profile.plot.log")
    script:
        "../scripts/plot_library_profile.R"