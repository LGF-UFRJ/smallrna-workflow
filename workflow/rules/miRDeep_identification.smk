
rule get_unassigned:
    input:
        lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.nh.bam.featureCounts.bam")
    output:
        os.path.join(mirdeep_dir, "{sample}.nh.unassigned.bam")
    log:
        os.path.join(mirdeep_dir, "log", "{sample}.nh.unassigned.log")
    conda:
        "../envs/pysam0.yaml"
    script:
        "../scripts/get_unassigned.py"

rule bam_to_collapsed_fasta:
    input:
        lambda wildcards: os.path.join(mirdeep_dir, wildcards.sample + ".nh.unassigned.bam")
    output:
        os.path.join(mirdeep_dir, "{sample}.reads_collapsed.fa")
    log:
        os.path.join(mirdeep_dir, "log", "{sample}.reads_collapsed.fa.log")
    conda:
        "../envs/mirdeep.yaml"
    shell:
        "samtools fasta {input} | python3 scripts/dedup_fasta.py | collapse_reads_md.pl - rpr > {output}"

rule map_mirdeep:
    input:
        lambda wildcards: os.path.join(mirdeep_dir, wildcards.sample + ".reads_collapsed.fa")
    output:
        os.path.join(mirdeep_dir, "{sample}.reads_collapsed.arf")
    log:
        os.path.join(mirdeep_dir, "log", "{sample}.reads_collapsed.arf.log")
    conda:
        "../envs/mirdeep.yaml"
    params:
        config["bowtie_vb_index"]
    shell:
        "mapper.pl {input} -c -p {params} -t {output}"

rule get_fasta_pre_miRNA:
    input:
        annot = config["annotations"],
        fa = config["genome_vb"]
    output:
        os.path.join(mirdeep_dir, "known_pre_miRNA.fa")
    log:
        os.path.join(mirdeep_dir, "log", "get_fasta_pre_miRNA.log")
    script:
        "../scripts/get_fasta_pre_miRNA.py"

rule format_genome:
    input:
        config["genome_vb"]
    output:
        os.path.join(mirdeep_dir, "genome.fa")
    script:
        "../scripts/format_genome.py"


rule miRDeep2:
    input:
        fasta = lambda wildcards: os.path.join(mirdeep_dir, wildcards.sample + ".reads_collapsed.fa"),
        arf = os.path.join(mirdeep_dir, "{sample}.reads_collapsed.arf"),
        genome = os.path.join(mirdeep_dir, "genome.fa"),
        known = os.path.join(mirdeep_dir, "known_pre_miRNA.fa"),
    log:
        # os.path.join(mirdeep_dir, "{sample}_miRDeep", "results.log")
        os.path.join(mirdeep_dir, "{sample}_miRDeep", "results.log")
    conda:
        "../envs/mirdeep.yaml"
    params:
        related_mature = config["related_mirs"],
        outdir = os.path.join(mirdeep_dir, "{sample}_miRDeep")
    shell:
        "mkdir -p {params.outdir} && "
        "cd {params.outdir} && "
        "miRDeep2.pl {input.fasta} {input.genome} "
        "{input.arf} none {params.related_mature} {input.known} 2> {log}"
