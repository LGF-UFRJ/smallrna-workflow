rule multiQC:
    output:
        report = os.path.join(multiqc_dir, "multiqc_report.html")
    params:
        f"-o {multiqc_dir}"
    shell:
        "multiqc {params} ../"