genome_size_file = os.path.join(bamcoverage_dir, "effective_genome_size.txt")
sfactors_file = os.path.join(bamcoverage_dir, "scallingFactors.tsv")

rule get_normalization_factors:
    input:
        expand(os.path.join(map_out_dir, "{sample}.sorted.counts.tsv"), sample = samplesheet["name"])
    output:
        os.path.join(bamcoverage_dir, "scallingFactors.tsv")
    shell:
        "python3 scripts/get_scale_factors.py -c {input} -o {output}"

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
        sfactors = os.path.join(bamcoverage_dir, "scallingFactors.tsv"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt"),
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.bw")
    params:
        main = f"-p 6 --effectiveGenomeSize $(cat {genome_size_file}) -bs 10",
        sf = sfactors_file,
        sample_name = lambda wildcards: wildcards.sample
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    threads: 6
    shell:
        "bamCoverage -b {input.bam} -o {output} {params.main} "
        "--scaleFactor $(cat {params.sf} | grep '{params.sample_name}' | cut -f 2)"

rule bamCoverage_uniquely_pos:
    input:
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        bai = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam.bai"),
        sfactors = os.path.join(bamcoverage_dir, "scallingFactors.tsv"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt")
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.pos.bw")
    params:
        main = f"-p 6 --effectiveGenomeSize $(cat {genome_size_file}) --filterRNAstrand reverse -bs 18",
        sf = sfactors_file,
        sample_name = lambda wildcards: wildcards.sample
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    threads: 6
    shell:
        "bamCoverage -b {input.bam} -o {output} {params.main} "
        "--scaleFactor $(cat {params.sf} | grep '{params.sample_name}' | cut -f 2)"

rule bamCoverage_uniquely_neg:
    input:
        bam = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam"),
        bai = lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".sorted.uniquely.bam.bai"),
        sfactors = os.path.join(bamcoverage_dir, "scallingFactors.tsv"),
        genome_size = os.path.join(bamcoverage_dir, "effective_genome_size.txt"),
    output:
        os.path.join(bamcoverage_dir, "{sample}.sorted.uniquely.neg.bw")
    params:
        main = f"-p 6 --effectiveGenomeSize $(cat {genome_size_file}) --filterRNAstrand forward -bs 18",
        sf = sfactors_file,
        sample_name = lambda wildcards: wildcards.sample
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    threads: 6
    shell:
        "bamCoverage -b {input.bam} -o {output} {params.main} "
        "--scaleFactor $(cat {params.sf} | grep '{params.sample_name}' | cut -f 2)"