te_dir = os.path.join(cluster_dir, "TEs/")

rule merge_clusters:
    input:
        expand(os.path.join(cluster_dir, "{sample}.clusters.nomir.bed"), sample = samplesheet["name"]),
    output:
        os.path.join(cluster_dir, "clusters.bed"),
    script:
        "../scripts/merge_clusters.py"

rule merge_stage_clusters:
    input:
        expand(os.path.join(cluster_dir, "{sample}.clusters.nomir.bed"), sample = samplesheet["name"]),
    output:
        expand(os.path.join(cluster_dir, "{stage}.clusters.nomir.bed"), stage = stages)
    script:
        "../scripts/merge_stage_clusters.py"

rule clusters_total:   
    input:
        bamfiles = expand(os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam"), sample = samplesheet["name"]),
        clusters = os.path.join(cluster_dir, "clusters.bed")
    output:
        summary = os.path.join(cluster_dir, "clusters.total.results.npz"),
        sfactors = os.path.join(cluster_dir, "clusters.total.scallingFactors.tsv"),
        outcount = os.path.join(cluster_dir, "clusters.counts.tsv"),
    threads: 6
    log:
        os.path.join(cluster_dir, "log", "clusters.multiBamSummary.log")
    shell:
        "multiBamSummary BED-file --BED {input.clusters} -p {threads} --outRawCounts {output.outcount} --scalingFactors {output.sfactors} --bamfiles {input.bamfiles} -o {output.summary}"

# rule clusters_uniquely:   
#     input:
#         bamfiles = expand(os.path.join(map_out_dir, "{sample}.sorted.uniquely.bam"), sample = samplesheet["name"]),
#         clusters = os.path.join(cluster_dir, "clusters.bed")
#     output:
#         summary = os.path.join(cluster_dir, "clusters.total.results.npz"),
#         sfactors = os.path.join(cluster_dir, "clusters.total.scallingFactors.tsv"),
#         outcount = os.path.join(cluster_dir, "clusters.counts.tsv"),
#     threads: 6
#     log:
#         os.path.join(cluster_dir, "log", "clusters.multiBamSummary.log")
#     shell:
#         "multiBamSummary BED-file --BED {input.clusters} -p {threads} --outRawCounts {output.outcount} --scalingFactors {output.sfactors} --bamfiles {input.bamfiles} -o {output.summary}"


rule get_pos_neg_piRNAs:   
    input:
        bam = lambda wildcards: os.path.join(cluster_pirna_dir, wildcards.sample + ".piRNAs.uniquely.bam"),
    output:
        pos = os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.pos.bam"),
        pos_bai = os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.pos.bam.bai"),
        neg = os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.neg.bam"),
        neg_bai = os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.neg.bam.bai"),
    threads: 6
    log:
        os.path.join(cluster_pirna_dir, "log", "{sample}.clusters.multiBamSummary.log")
    shell:
        "samtools view -b -F 16 {input.bam} > {output.pos} && "
        "samtools view -b -f 16 {input.bam} > {output.neg} && "
        "samtools index {output.pos} && " 
        "samtools index {output.neg}"

rule clusters_multiBamSummary:   
    input:
        pos = expand(os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.pos.bam"), sample = samplesheet["name"]),
        neg = expand(os.path.join(cluster_pirna_dir, "{sample}.piRNAs.uniquely.neg.bam"), sample = samplesheet["name"]),
        clusters = os.path.join(cluster_dir, "clusters.bed"),
    output:
        pos = os.path.join(cluster_dir, "clusters.pos.counts.tsv"),
        neg = os.path.join(cluster_dir, "clusters.neg.counts.tsv"),
        pos_npz = os.path.join(cluster_dir, "clusters.pos.counts.npz"),
        neg_npz = os.path.join(cluster_dir, "clusters.neg.counts.npz"),
    threads: 6
    log:
        os.path.join(cluster_dir, "log", "clusters.multiBamSummary.log")
    shell:
        "multiBamSummary BED-file --BED {input.clusters} -p 6 --outRawCounts {output.pos} --bamfiles {input.pos} -o {output.pos_npz} && "
        "multiBamSummary BED-file --BED {input.clusters} -p 6 --outRawCounts {output.neg} --bamfiles {input.neg} -o {output.neg_npz}"

rule clusters_info:   
    input:
        os.path.join(cluster_dir, "clusters.counts.tsv"),
        os.path.join(cluster_dir, "clusters.pos.counts.tsv"),
        os.path.join(cluster_dir, "clusters.neg.counts.tsv"),
        expand(os.path.join(map_out_dir, "{sample}.sorted.counts.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(cluster_dir, "clusters.final.tsv"),
    threads: 6
    log:
        os.path.join(cluster_dir, "log", "clusters.multiBamSummary.log")
    script:
        "../scripts/clusters_rpm.py"

rule final_to_bed:
    input:
        os.path.join(cluster_dir, "clusters.final.tsv")
    output:
        os.path.join(cluster_dir, "clusters.final.bed")
    shell:
        "bash scripts/cluster_fmt_final.sh {input} > {output}"

rule bedToBigBed:
    input:
        clusters = os.path.join(cluster_dir, "clusters.final.tsv"),
        genome_index = config["genome_index"]
    output:
        sortedcl = os.path.join(cluster_dir, "clusters.sorted.bed"),
        bb = os.path.join(cluster_dir, "clusters.bb"),
    shell:
        "bash scripts/clusters_bigBed.sh {input.clusters} {input.genome_index} {output.sortedcl} {output.bb}"

rule clusters_overlap:
    input:
        clusters = os.path.join(cluster_dir, "clusters.bed"),
        files = expand(os.path.join(cluster_dir, "{sample}.clusters.nomir.bed"), sample = samplesheet["name"])
    output:
        os.path.join(cluster_dir, "clusters.overlap.tsv")
    shell:
        "bash scripts/clusters_overlap.sh {input.clusters} {input.files} >> {output}"

#rule overlap_TE_top_clusters:
    #input:
        #clusters = os.path.join(cluster_dir, "top10_clusters.bed"),
        #tes = os.path.join(cluster_pirna_dir, "TE_annotation.bed")
    #output:
        #os.path.join(te_dir, "top10_clusters.te.overlap.bed")
    #shell:
        #"bedtools intersect -wo -a {input.top_cl} -b {input.tes} > {output}"

rule clusters_top10_bed:
    input:
        os.path.join(cluster_dir, "clusters.final.tsv"),
    output:
        os.path.join(cluster_dir, "top10_clusters.bed")
    shell:
        "bash scripts/get_top10_clusters.sh {input} {output}"

rule overlap_TE_top_clusters:
    input:
        top_cl = os.path.join(cluster_dir, "top10_clusters.bed"),
        tes = os.path.join(cluster_pirna_dir, "TE_annotation.bed")
    output:
        os.path.join(te_dir, "top10_clusters.te.overlap.bed")
    shell:
        "bedtools intersect -wo -a {input.top_cl} -b {input.tes} > {output}"

rule cluster_TE_count:
    input:
        os.path.join(te_dir, "top10_clusters.te.overlap.bed")
    params:
        te_dir
    output:
        [os.path.join(te_dir, f"RPCL{i}.TEs.count.tsv") for i in range(1,11)]
    shell:
        "bash scripts/get_count_per_cluster.sh {input} {params}"

rule clusters_TE_strand_count:
    input:
        sortedcl = os.path.join(cluster_dir, "clusters.sorted.bed"),
        tes = os.path.join(cluster_pirna_dir, "TE_annotation.bed")
    output:
        os.path.join(cluster_dir, "clusters.te.strand.count.tsv"),
    shell:
        "bash scripts/clusters_TE_strand_count.sh {input.tes} {input.sortedcl} > {output}"

rule bam_to_bed_uniquely:
    input: 
        lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam")
    output: 
        os.path.join(map_out_dir, "{sample}.sorted.uniquely.bed")
    shell: 
        "bedtools bamtobed -i {input} > {output}"

rule merge_uniquely_bed:
    input: 
        expand(os.path.join(map_out_dir, "{sample}.sorted.uniquely.bed"), sample = samplesheet["name"])
    output: 
        os.path.join(map_out_dir, "all.sorted.uniquely.bed")
    shell: 
        "cat {input} | bedtools sort > {output}"

rule top10_pos_coverage:
    input: 
        all_bed = os.path.join(map_out_dir, "all.sorted.uniquely.bed"),
        top10 = os.path.join(cluster_dir, "top10_clusters.bed")
    output: 
        os.path.join(cluster_dir, "top10_clusters.cov.pos.tsv")
    shell: 
        "bedtools coverage -d -s -a <(sed 's/\./\+/' {input.top10}) -b {input.all_bed} > {output}"

rule top10_neg_coverage:
    input: 
        all_bed = os.path.join(map_out_dir, "all.sorted.uniquely.bed"),
        top10 = os.path.join(cluster_dir, "top10_clusters.bed")
    output: 
        os.path.join(cluster_dir, "top10_clusters.cov.neg.tsv")
    shell: 
        "bedtools coverage -d -S -a <(sed 's/\./\+/' {input.top10}) -b {input.all_bed} > {output}"

rule format_bed:
    input: 
        os.path.join(cluster_dir, "clusters.final.tsv"),
    output: 
        os.path.join(cluster_dir, "clusters.fmt.bed"),
    shell: 
        "bash scripts/cl_cfmt.sh {input} > {output}"