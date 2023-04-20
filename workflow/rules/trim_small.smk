rule fastq_files:
    input:
        config["samplesheet"]
    output:
        expand(os.path.join(fastq_dir, "{sample}.fastq.gz"), sample = samplesheet["name"])
    params:
        fastq_dir
    log:
        os.path.join(fastq_dir, "log", "fastq_files.log")
    script:
        "../scripts/get_input_fastqs.py"

rule trim_small:
    input:
        expand(os.path.join(fastq_dir, "{sample}.fastq.gz"),
               sample = samplesheet["name"])
    output:
        expand(os.path.join(trim_outdir, "{sample}_trimmed.fq.gz"), sample = samplesheet["name"]),
    priority: 100
    log:
        os.path.join(trim_outdir, "log", "trim_small.log")
    params:
        f"--fastqc --length 18 --max_length 35 -o {trim_outdir}"
    shell:
        "trim_galore {params} {input}"

