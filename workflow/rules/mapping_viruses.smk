
index_vir = config["bowtie_virus_index"]


rule bowtie_viruses:
    input:
        lambda wildcards: os.path.join(trim_outdir, wildcards.sample + "_trimmed.fq.gz")
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


rule length_distribution_viruses:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.lengths.tsv")
    params:
        os.path.join(map_out_dir, "{sample}.sorted.counts.tsv")
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

# -----------------------------------------------------------------------------


rule sep_sense_antisense_viruses:
    input:
        bam = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
    output:
        sense = os.path.join(map_out_virus_dir, "{sample}.viruses.sense.bam"),
        antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.antisense.bam"),
    shell:
        "samtools view -Sb -F 16 {input.bam} > {output.sense} && "
        "samtools view -Sb -f 16 {input.bam} > {output.antisense}"


rule bamtobed_sense_viruses:
    input:
        bam = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sense.bam"),
        #antisense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.antisense.bam"),
    output:
        bed = os.path.join(map_out_virus_dir, "{sample}.viruses.sense.bed"),
        #antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.antisense.bed"),
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed}"

rule bamtobed_antisense_viruses:
    input:
        #lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sense.bam"),
        bam = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.antisense.bam"),
    output:
        bed = os.path.join(map_out_virus_dir, "{sample}.viruses.antisense.bed"),
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed}"

rule count_overlaps_viruses:
    input:
        sense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sense.bed"),
        antisense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.antisense.bed"),
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.gz"),
    shell:
        "bedtools intersect -wo -S -a {input.sense} -b {input.antisense} | gzip - > {output}"


rule filter_5p_overlap_viruses:
    input:
        os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.gz"),
    output:
        os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.gz"),
    shell:
        "bash scripts/get_5p_overlap.sh {input} > {output}"

rule separate_5p_bed_sense_viruses:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.gz"),
    output:
        sense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.bed"),
    shell:
        "bash scripts/get_5p_bed.sh {input} sense > {output.sense}"

rule separate_5p_bed_antisense_viruses:
    input:
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.gz"),
    output:
        antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.bed"),
    shell:
        "bash scripts/get_5p_bed.sh {input} antisense > {output.antisense}"

rule get_5p_fasta_viruses:
    input:
        genome = config["genome_virus"],
        sense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.sense.bed"),
        antisense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.antisense.bed"),
    output:
        sense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.fa"),
        antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.fa"),
    shell:
        "bedtools getfasta -s -name -fi {input.genome} -bed {input.sense} > {output.sense} && "
        "bedtools getfasta -s -name -fi {input.genome} -bed {input.antisense} > {output.antisense}"

rule get_seqs_logo_viruses:
    input:
        sense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.sense.fa"),
        antisense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.antisense.fa"),
    output:
        sense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.seqs.txt"),
        antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.seqs.txt"),
    shell:
        "bash scripts/get_seqs_for_logo.sh {input.sense} > {output.sense} && "
        "bash scripts/get_seqs_for_logo.sh {input.antisense} > {output.antisense}"

rule merge_seqs_viruses:
    input: 
        sense = expand(os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.seqs.txt"), sample = samplesheet["name"]),
        antisense = expand(os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.seqs.txt"), sample = samplesheet["name"]),
    output: 
        sense = os.path.join(map_out_virus_dir, "all.sense.seqs.txt"),
        antisense = os.path.join(map_out_virus_dir, "all.antisense.seqs.txt"),
    shell: 
        "cat {input.sense} > {output.sense} && "
        "cat {input.antisense} > {output.antisense}"

rule plot_seq_logo_viruses:
    input:
        # sense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.sense.seqs.txt"),
        # antisense = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.5p.antisense.seqs.txt"),
        sense = os.path.join(map_out_virus_dir, "all.sense.seqs.txt"),
        antisense = os.path.join(map_out_virus_dir, "all.antisense.seqs.txt"),
        # sense = expand(os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.seqs.txt"), sample = samplesheet["name"]),
        # antisense = expand(os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.seqs.txt"), sample = samplesheet["name"]),
    output:
        # sense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.sense.seqs.logo.png"),
        # antisense = os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.5p.antisense.seqs.png"),
        sense = os.path.join(map_out_virus_dir, "all.viruses.overlap.tsv.5p.sense.seqs.logo.png"),
        antisense = os.path.join(map_out_virus_dir, "all.viruses.overlap.tsv.5p.antisense.seqs.png"),

    shell:
        "Rscript scripts/plot_seqlogo.R {input.sense} {output.sense} && "
        "Rscript scripts/plot_seqlogo.R {input.antisense} {output.antisense}"

rule nolowcomp_viruses:
    output:
        os.path.join(map_out_virus_dir, "nolowcomp")
    shell:
        "touch {output}"

rule merge_overlaps_viruses:
    input:
        expand(os.path.join(map_out_virus_dir, "{sample}.viruses.overlap.tsv.gz"), sample = samplesheet["name"]),
    output:
        os.path.join(map_out_virus_dir, "all.viruses.overlap.tsv.gz")
    shell:
        "cat {input} > {output}"

rule merge_mrc_viruses:
    input:
        expand(os.path.join(map_out_vb_dir, "{sample}.mapreadcount.tsv"), sample = samplesheet["name"]),
    output:
        os.path.join(map_out_virus_dir, "all.mapreadcount.tsv")
    shell:
        "cat {input} > {output}"

rule pingpong_signal_viruses:
    input:
        lowcomp = "nolowcomp",
        # mrc = lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".mapreadcount.tsv"),
        # overlap = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.overlap.tsv.gz"),
        mrc = os.path.join(map_out_virus_dir, "all.mapreadcount.tsv"),
        overlap = os.path.join(map_out_virus_dir, "all.viruses.overlap.tsv.gz"),
    output:
        # os.path.join(map_out_virus_dir, "{sample}.viruses.pingpong.tsv"),
        os.path.join(map_out_virus_dir, "all.viruses.pingpong.tsv"),
    shell:
        "python3 scripts/pingpong_signal.py {input.lowcomp} {input.mrc} {input.overlap} > {output}"

rule count_viruses:
    input: 
        viruses = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
        mrc = lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".mapreadcount.tsv"),
    output: 
        os.path.join(map_out_virus_dir, "{sample}.virus-count.tsv")
    shell:
        "python3 scripts/count.py {input.viruses} {input.mrc} > {output}"



rule coverage_pos_viruses:
    input: 
        viruses = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
    output: 
        os.path.join(map_out_virus_dir, "{sample}.genomecov.rpm.pos.tsv")
    params:
        sf = sfactors_file,
        sample_name = lambda wildcards: wildcards.sample
    shell: 
        "bedtools genomecov -d -strand + -scale $(cat {params.sf} | grep '{params.sample_name}' | cut -f 2) "
        "-ibam {input.viruses} > {output}"

rule coverage_neg_viruses:
    input: 
        viruses = lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
    output: 
        os.path.join(map_out_virus_dir, "{sample}.genomecov.rpm.neg.tsv")
    params:
        sf = sfactors_file,
        sample_name = lambda wildcards: wildcards.sample
    shell: 
        "bedtools genomecov -d -strand '-' -scale $(cat {params.sf} | grep '{params.sample_name}' | cut -f 2) "
        "-ibam {input.viruses} > {output}"

# -----------------------------------------------------------------------------

rule lendist_each_virus:
    input: 
        lambda wildcards: os.path.join(map_out_virus_dir, wildcards.sample + ".viruses.sorted.bam"),
    output: 
        os.path.join(map_out_virus_dir, "{sample}.viruses.sorted.lendist.each.tsv")
    shell: 
        "python3 scripts/lendist_each_viruses.py {input} {output}"