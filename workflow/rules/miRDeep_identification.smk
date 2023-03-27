
rule get_repeats_unassigned:
    input:
        lambda wildcards: os.path.join(pirna_dir, wildcards.sample + ".genes.unassigned.bam.featureCounts.bam")
    output:
        os.path.join(unassigned_dir, "{sample}.gene.repeats.unassigned.bam")
    log:
        os.path.join(unassigned_dir, "log", "{sample}.genome.repeats.unassigned.log")
    conda:
        "../envs/pysam0.yaml"
    script:
        "../scripts/get_unassigned.py"

rule bam_to_collapsed_fasta:
    input:
        lambda wildcards: os.path.join(unassigned_dir, wildcards.sample + ".gene.repeats.unassigned.bam")
    output:
        os.path.join(unassigned_dir, "{sample}.reads_collapsed.fa")
    log:
        os.path.join(unassigned_dir, "log", "{sample}.reads_collapsed.fa.log")
    conda:
        "../envs/mirdeep.yaml"
    shell:
        "samtools fasta {input} | python3 scripts/dedup_fasta.py | collapse_reads_md.pl - rpr > {output}"

rule map_mirdeep:
    input:
        lambda wildcards: os.path.join(unassigned_dir, wildcards.sample + ".reads_collapsed.fa")
    output:
        os.path.join(unassigned_dir, "{sample}.reads_collapsed.arf")
    log:
        os.path.join(unassigned_dir, "log", "{sample}.reads_collapsed.arf.log")
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
        os.path.join(unassigned_dir, "known_pre_miRNA.fa")
    log:
        os.path.join(unassigned_dir, "log", "get_fasta_pre_miRNA.log")
    script:
        "../scripts/get_fasta_pre_miRNA.py"

rule format_genome:
    input:
        config["genome_vb"]
    output:
        os.path.join(unassigned_dir, "genome.fa")
    script:
        "../scripts/format_genome.py"


rule miRDeep2:
    input:
        lambda wildcards: os.path.join(unassigned_dir, wildcards.sample + ".reads_collapsed.fa")
    output:
        # os.path.join(unassigned_dir, "{sample}_miRDeep", "results.log")
        os.path.join(unassigned_dir, "{sample}_miRDeep", "done.txt")
    conda:
        "../envs/mirdeep.yaml"
    params:
        genome = os.path.join(unassigned_dir, "genome.fa"),
        arf = os.path.join(unassigned_dir, "{sample}.reads_collapsed.arf"),
        related_mature = config["related_mirs"],
        known = os.path.join(unassigned_dir, "known_pre_miRNA.fa"),
        outdir = os.path.join(unassigned_dir, "{sample}_mirDeep")
    shell:
        "mkdir -p {params.outdir} && "
        "cd {params.outdir} && "
        "miRDeep2.pl {input} {params.genome} "
        "{params.arf} none {params.related_mature} {params.known} 2> results.log && "
        "touch done.txt"
