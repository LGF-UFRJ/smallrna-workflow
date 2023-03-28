rule link_inputs:
    input:
        config["samplesheet"]
    output:
        expand(os.path.join(trim_outdir, "file_links", "{sample}.fastq.gz"), sample = samplesheet["name"])
    params:
        os.path.join(trim_outdir, "file_links/")
    log:
        os.path.join(trim_outdir, "log", "link_inputs.log")
    script:
        "../scripts/get_input_links.py"

rule trim_small:
    input:
        expand(os.path.join(trim_outdir, "file_links", "{sample}.fastq.gz"),
               sample = samplesheet["name"])
    output:
        os.path.join(trim_outdir, "{sample}_trimmed.fq.gz")
    log:
        os.path.join(trim_outdir, "log", "{sample}.trim_small.log")
    params:
        f"--fastqc --length 18 --max_length 35 -o {trim_outdir}"
    shell:
        "trim_galore {params} {input}"

