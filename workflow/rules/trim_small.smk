



# rule link_inputs:
#     input:
#         # expand(os.path.join(basepath, "{sample}.fastq.gz"), sample = basenames_n_dir)
#         lambda wildcards: os.path.join(map_out_dir, wildcards.sample + ".fastq.gz")
#     output:
#         outfiles = expand(os.path.join(trim_outdir, "file_links", "{filename}.fastq.gz"), filename = basenames)
#     params:
#         outdir = os.path.join(trim_outdir, "file_links/")
#     log:
#         os.path.join(trim_outdir, "log", "link_inputs.log")
#     shell:
#         "ln -s {input} {params.outdir}"



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
        expand(os.path.join(trim_outdir, "file_links", "{filename}.fastq.gz"),
               filename = basenames)
    output:
        os.path.join(trim_outdir, "{filename}_trimmed.fq.gz")
    log:
        os.path.join(trim_outdir, "log", "{filename}.trim_small.log")
    params:
        f"--fastqc --length 18 --max_length 35 -o {trim_outdir}"
    shell:
        "trim_galore {params} {input}"

