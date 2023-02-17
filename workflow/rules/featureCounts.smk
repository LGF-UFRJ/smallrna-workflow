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
        os.path.join(featureCount_dir, "featureCounts.counts.tsv")
    log:
        os.path.join(featureCount_dir, "log", "featureCounts.log")
    threads: 6
    params:
        " -F SAF -M -O -R CORE --fraction"
    shell:
        "featureCounts {params} -T {threads} -a {input.saf} -o {output} {input.bam}"


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