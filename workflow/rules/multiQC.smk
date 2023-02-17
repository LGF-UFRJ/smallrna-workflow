rule multiQC:
    input:
        os.path.join(deeptools_dir, "plotCorrelation.png")
    output:
        os.path.join(multiqc_dir, "multiqc_report.html")
    params:
        f"-o {multiqc_dir}"
    shell:
        "multiqc {params} ../"