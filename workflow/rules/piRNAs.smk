rule get_piRNAs:
    input:
        lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.nh.bam.featureCounts.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.bam")
    script:
        "../scripts/get_piRNAs.py"

rule pingpong:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.pingpong.tsv"),
        os.path.join(pirna_dir, "{sample}.piRNAs.pairs.tsv"),
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

rule split_piRNA_pairs:
    input:
        pairs = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.pairs.tsv"),
        sense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.bam"),
        antisense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.bam"),
    output:
        sense = os.path.join(pirna_dir, "{sample}.piRNAs.sense.pairs.bam"),
        antisense = os.path.join(pirna_dir, "{sample}.piRNAs.antisense.pairs.bam"),
    shell:
        "bash scripts/split_piRNAs.sh {input.pairs} {input.sense} {input.antisense} {output.sense} {output.antisense} &> out.log"

rule get_sense_fa:
    input:
        sense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.pairs.bam"),
        antisense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.pairs.bam"),
    output:
        sense = os.path.join(pirna_dir, "{sample}.piRNAs.sense.pairs.fa"),
        antisense = os.path.join(pirna_dir, "{sample}.piRNAs.antisense.pairs.fa"),
    shell:
        "bash scripts/get_sense_fasta.sh {input.sense} > {output.sense} &&"
        "bash scripts/get_sense_fasta.sh {input.antisense} > {output.antisense}"

rule get_seqs_logo:
    input:
        sense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.pairs.fa"),
        antisense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.pairs.fa"),
    output:
        sense = os.path.join(pirna_dir, "{sample}.piRNAs.sense.pairs.seqs.txt"),
        antisense = os.path.join(pirna_dir, "{sample}.piRNAs.antisense.pairs.seqs.txt"),
    shell:
        "bash scripts/get_seqs_for_logo.sh {input.sense} > {output.sense} &&"
        "bash scripts/get_seqs_for_logo.sh {input.antisense} > {output.antisense}"

rule plot_seq_logo:
    input:
        sense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.pairs.seqs.txt"),
        antisense = lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.pairs.seqs.txt"),
    output:
        sense = os.path.join(pirna_dir, "{sample}.piRNAs.sense.pairs.seqs.logo.png"),
        antisense = os.path.join(pirna_dir, "{sample}.piRNAs.antisense.pairs.seqs.logo.png"),
    shell:
        "Rscript scripts/plot_seqlogo.R {input.sense} {output.sense} && "
        "Rscript scripts/plot_seqlogo.R {input.antisense} {output.antisense}"

rule nt_bias_sense:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.sense.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.sense.ntfreq.tsv")
    script:
        "../scripts/nucleotide_bias.py"


rule nt_bias_antisense:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".piRNAs.antisense.bam")
    output:
        os.path.join(pirna_dir, "{sample}.piRNAs.antisense.ntfreq.tsv")
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

#rule TE_piRNA_quantification:
    #input:
        #repeat_saf = os.path.join(featureCount_dir, "annotations", "repeats.saf"),
        #bam = expand(os.path.join(map_out_vb_dir, "{sample}.sorted.nh.bam"), sample = samplesheet["name"])
    #output:
        #counts = os.path.join(featureCount_dir, "featureCounts.counts.tsv"),
        #outbam = expand(os.path.join(featureCount_dir, "{sample}.sorted.nh.bam.featureCounts.bam"), sample = samplesheet["name"])
    #log:
        #os.path.join(featureCount_dir, "log", "featureCounts.log")
    #threads: 6
    #params:
        #" -F SAF -M -O -R BAM --fraction"
    #shell:
        #"featureCounts {params} -T {threads} -a {input.saf} -o {output.counts} {input.bam}"

rule get_mapreadcount_counts:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.bam")
    output:
        os.path.join(map_out_vb_dir, "{sample}.mapreadcount.tsv")
    shell:
        "bash scripts/mapcount.sh {input} > {output}"

rule get_VB_TEs_bed:
    input: 
        config["repeat_annotation"]
    output:
        os.path.join(pirna_dir, "VB_TE_annotation.bed")
    shell:
        "bash scripts/get_TE_bed.sh {input} > {output}"

rule count_te_piRNAs:
    input:
        mrc = lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".mapreadcount.tsv"),
        fc = lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.nh.bam.featureCounts.bam"),
    output:
        os.path.join(pirna_dir, "{sample}.te.piRNA.count.tsv"),
    shell:
        "bash scripts/count_te.sh {input.mrc} {input.fc} > {output}"

