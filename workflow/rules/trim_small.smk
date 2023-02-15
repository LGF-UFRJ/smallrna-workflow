

rule link_inputs:
    input:
        expand(os.path.join(basepath, "{sample}.fastq.gz"), sample = basenames_n_dir)
    output:
        outfiles = expand(os.path.join(trim_outdir, "file_links", "{filename}.fastq.gz"), filename = basenames)
    params:
        outdir = os.path.join(trim_outdir, "file_links/")
    shell:
        "ln -s {input} {params.outdir}"

rule trim_small:
    input:
        expand(os.path.join(trim_outdir, "file_links", "{filename}.fastq.gz"),
               filename = basenames)
    output:
        os.path.join(trim_outdir, "{filename}_trimmed.fq.gz")
    log:
        os.path.join(trim_outdir, "log", "{filename}.log.txt")
    params:
        f"--fastqc --length 18 --max_length 35 -o {trim_outdir}"
    shell:
        "trim_galore {params} {input} &> {log}"

