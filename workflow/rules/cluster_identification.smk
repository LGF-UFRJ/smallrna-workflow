hic_fa = config["genome"]
rm_out = config["repeat_masker_out"]


rule bam2sam:
    input:
        lambda wildcards: os.path.join(map_out_vb_dir, wildcards.sample + ".sorted.uniquely.bam")
    output:
        os.path.join(cluster_dir, "{sample}.sorted.uniquely.sam")
    shell:
        "samtools view -h {input} > {output}"


rule proTRAC:
    input:
        lambda wildcards: os.path.join(cluster_dir, wildcards.sample + ".sorted.uniquely.sam")
    log:
        os.path.join(cluster_dir, "{sample}_clusters", "protrac_results.log")
    params:
        outdir = os.path.join(cluster_dir, "{sample}_clusters"),
        main = f"-g {hic_fa} -repeatmasker {rm_out} -format SAM -1Tor10A 0.0 -distr 1-100"
    conda:
        "../envs/protrac.yaml"
    shell:
        "cd {params.outdir} && "
        "proTRAC_2.4.2.pl -m {input} {params.main} 2> {log}"
