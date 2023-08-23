rule get_piRNAs:
    input:
        lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.nh.bam.featureCounts.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.bam")
    #conda:
        #"../envs/pysam0.yaml"
    script:
        "../scripts/get_piRNAs.py"

rule pingpong:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.pingpong.tsv"),
        os.path.join(pirna_dir, "{sample}.piRNAs.pairs.tsv"),
    #conda:
        #"../envs/pysam0.yaml"
    script:
        "../scripts/pingpong.py"

rule plot_pingpong:
    input:
        expand(os.path.join(pirna_dir, "{sample}.piRNAs.pingpong.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(pirna_dir, "ridges.piRNAs.pingpong.plot.png"),
        os.path.join(pirna_dir, "lines.piRNAs.pingpong.plot.png"),
    script:
        "../scripts/plot_pingpong.R"

rule sep_sense_antisense:
    input:
        pirnas = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.bam"),
        gtf = config["repeat_annotation"]
    output:
        sense = os.path.join(pirna_dir, "{sample}.piRNAs.sense.bam"),
        antisense = os.path.join(pirna_dir, "{sample}.piRNAs.antisense.bam"),
    shell:
        "samtools view -b {input.pirnas} | bedtools intersect -s -abam stdin -b {input.gtf} > {output.sense} && "
        "samtools view -b {input.pirnas} | bedtools intersect -S -abam stdin -b {input.gtf} > {output.antisense}"


rule nt_bias_sense:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.sense.ntfreq.tsv")
    #conda:
        #"../envs/pysam0.yaml"
    script:
        "../scripts/nucleotide_bias.py"


rule nt_bias_antisense:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.antisense.ntfreq.tsv")
    #conda:
        #"../envs/pysam0.yaml"
    script:
        "../scripts/nucleotide_bias.py"

rule plot_nt_freq_sense_piRNA:
    input:
        expand(os.path.join(pirna_dir, "{sample}.piRNAs.sense.ntfreq.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(pirna_dir, "piRNAs.sense.ntfreq.plot.png"),
    log:
        os.path.join(pirna_dir, "log", "sense.ntfreq.plot.log")
    script:
        "../scripts/plot_ntfreq.R"

rule plot_nt_freq_antisense_piRNA:
    input:
        expand(os.path.join(pirna_dir, "{sample}.piRNAs.antisense.ntfreq.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(pirna_dir, "piRNAs.antisense.ntfreq.plot.png"),
    log:
        os.path.join(pirna_dir, "log", "antisense.ntfreq.plot.log")
    script:
        "../scripts/plot_ntfreq.R"


rule length_distribution_piRNAs:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.lengths.tsv")
    params:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
    #conda:
        #"../envs/pysam0.yaml"
    log:
        os.path.join(pirna_dir, "log", "{sample}.length_distribution_piRNAs.log")
    script:
        "../scripts/length_distribution.py"

rule plot_length_distribution_piRNA:
    input:
        expand(os.path.join(pirna_dir, "{sample}.piRNAs.lengths.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(pirna_dir, "piRNAs.ridges.lengths.plot.png"),
        os.path.join(pirna_dir, "piRNAs.line.lengths.plot.png"),
    log:
        os.path.join(pirna_dir, "log", "lengths.plot.log")
    script:
        "../scripts/plot_length_distribution.R"

rule TE_piRNA_quantification:
    input:
        repeat_saf = os.path.join(featureCount_dir, "annotations", "repeats.saf"),
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


