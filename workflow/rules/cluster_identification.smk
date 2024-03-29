hic_fa = config["genome"]
rm_out = config["repeat_masker_out"]
cluster_pirna_dir = os.path.join(cluster_dir, "piRNAs")

localrules: get_uniquely_piRNAs

rule get_TEs_bed:
    input: 
        config["repeats_hic"]
    output:
        os.path.join(cluster_pirna_dir, "TE_annotation.bed")
    shell:
        "bash scripts/get_TE_bed.sh {input} > {output}"

rule get_uniquely_piRNAs:
    input: 
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        tes = os.path.join(cluster_pirna_dir, "TE_annotation.bed")
    output: 
        os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.bam")
    shell:
        "bash scripts/get_hic_piRNAs.sh {input.bam} {input.tes} {output}"

# rule get_piRNA_ids:
#     input:
#         lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.bam"),
#     output:
#         os.path.join(cluster_pirna_dir, "{sample}.uniquely.ids.txt")
#     log:
#         os.path.join(cluster_pirna_dir, "log", "{sample}.get_piRNA_ids.log")
#     shell:
#         "bash scripts/get_unique_ids.sh {input} {output}"

rule get_sliding_window:
    input:
        config["genome_index"]
    output:
        os.path.join(cluster_dir, "windows.bed")
    shell:
        "bedtools makewindows -w 5000 -s 1000 -g {input} > {output}"


# rule get_uniquely_piRNAs:
#     input:
#         bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.bam"),
#         ids = lambda wildcards: os.path.join(cluster_pirna_dir, wildcards.sample + ".uniquely.ids.txt")
#     params:
#         sam = lambda wildcards: os.path.join(cluster_pirna_dir, wildcards.sample + ".piRNAs.sam"),
#     output:
#         os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.bam"),
#     log:
#         os.path.join(cluster_pirna_dir, "log", "{sample}.get_uniquely_piRNAs.log")
#     shell:
#         "samtools view -H {input.bam} > {params.sam} && "
#         "samtools view {input.bam} | grep -Ff {input.ids} >> {params.sam} && "
#         "samtools view -Sb {input.bam} > {output} && "
#         "samtools index {output} && "
#         "rm {params.sam}"



rule windows_multiBamSummary:   
    input:
        expand(os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.bam"), sample = samplesheet["name"]),
        # expand(os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam"), sample = samplesheet["name"]),
    output:
        summary = os.path.join(cluster_dir, "windows.results.npz"),
        outcount = os.path.join(cluster_dir, "windows.counts.tsv")
    threads: 6
    shell:
        "multiBamSummary bins -bs 1000 -p {threads} --outRawCounts {output.outcount} --bamfiles {input} -o {output.summary}"

rule normalize_windows:
    input:
        windows = os.path.join(cluster_dir, "windows.counts.tsv"),
        count_files = expand(os.path.join(map_out_dir, "{sample}.sorted.counts.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(cluster_dir, "windows.norm.tsv")
    shell:
      "python3 scripts/windows_norm.py -w {input.windows} -c {input.count_files} -o {output}"



rule extract_clusters:
    input:
        bam = lambda wildcards: os.path.join(cluster_pirna_dir, wildcards.sample + ".piRNAs.uniquely.bam"),
        windows = os.path.join(cluster_dir, "windows.bed")
    output:
        os.path.join(cluster_dir, "{sample}.clusters.bed"), 
        os.path.join(cluster_dir, "{sample}.window_count.tsv"), 
        os.path.join(cluster_dir, "{sample}.merged.bed"), 
    shell:
        "bash scripts/find_clusters.sh {input.bam} {input.windows}"


# rule extract_clusters:
#     input:
#         os.path.join(cluster_dir, "windows.norm.tsv")
#     params:
#         cluster_dir
#     output:
#         expand(os.path.join(cluster_dir, "{sample}.clusters.bed"), sample = samplesheet["name"])
#     shell:
#         "bash scripts/clusters.sh {input} {params}"

rule vb_bigBedToBed:
    input:
        config["genome_vb_transfered"],
    output:
        os.path.join(cluster_dir, "VB_transfered_annotation.bed")
    conda:
        "../envs/ucsc_tools.yaml"
    shell:
        "/home5/attilio/programs/ucsc_tools/bigBedToBed {input} {output}"

rule get_miRNAs_ids:
    input:
        os.path.join(featureCount_dir, "annotations", "names.tsv"),
    output:
        os.path.join(cluster_dir, "miRNA_ids.txt")
    shell:
        "cat {input} | grep -i mir | cut -f 1 > {output}"

rule get_miRNAs_bed:
    input:
        annotation = os.path.join(cluster_dir, "VB_transfered_annotation.bed"),
        mir_ids = os.path.join(cluster_dir, "miRNA_ids.txt")
    output:
        os.path.join(cluster_dir, "miRNAs.bed")
    shell:
        "cat {input.annotation} | grep -Ff {input.mir_ids} > {output}"

rule filter_clusters:
    input:
        clusters = lambda wildcards: os.path.join(cluster_dir, wildcards.sample + ".clusters.bed"),
        mirnas = os.path.join(cluster_dir, "miRNAs.bed"),
    output:
        os.path.join(cluster_dir, "{sample}.clusters.nomir.bed"),
    shell:
        "bedtools intersect -v -a {input.clusters} -b {input.mirnas} > {output}"
