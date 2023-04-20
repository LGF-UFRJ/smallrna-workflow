genome_size_file = os.path.join(bamcoverage_dir, "effective_genome_size.txt")

rule get_effective_genome_size:
    input:
        config["genome"]
    output:
        os.path.join(bamcoverage_dir, "effective_genome_size.txt")
    shell:
        "faCount {input} | tail -n 1 | cut -f 2 > {output}"


rule bamCoverage_uniquely:
    input:
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        bai = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam.bai"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt")
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.bw")
    params:
        f"-p 6 --normalizeUsing CPM --effectiveGenomeSize $(cat {genome_size_file})"
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    shell:
        "bamCoverage -b {input.bam} -o {output} {params}"


rule bamCoverage_uniquely_pos:
    input:
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        bai = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam.bai"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt")
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.pos.bw")
    params:
        f"-p 6 --normalizeUsing CPM --effectiveGenomeSize $(cat {genome_size_file}) "
        "--filterRNAstrand reverse"
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    shell:
        "bamCoverage -b {input.bam} -o {output} {params}"

rule bamCoverage_uniquely_neg:
    input:
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        bai = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam.bai"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt"),
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.neg.bw")
    params:
        f"-p 6 --normalizeUsing CPM --effectiveGenomeSize $(cat {genome_size_file}) "
        "--filterRNAstrand forward -bs 20"
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    shell:
        "bamCoverage -b {input.bam} -o {output} {params}"
