
rule multiBamSummary:   
    input:
        expand(os.path.join(map_out_dir, "{sample}.sorted.bam"), sample = basenames)
    output:
        summary = os.path.join(deeptools_dir, "multiBamSummary.results.npz"),
        sfactors = os.path.join(deeptools_dir, "multiBamSummary.scallingFactors.tsv")
    threads: 6
    log:
        os.path.join(deeptools_dir, "log", "multiBamSummary.log")
    shell:
        "multiBamSummary bins -p {threads} --scalingFactors {output.sfactors} --bamfiles {input} -o {output.summary}"


rule plotCorrelation:
    input:
        os.path.join(deeptools_dir, "multiBamSummary.results.npz")
    output:
        corr = os.path.join(deeptools_dir, "plotCorrelation.png"),
        table = os.path.join(deeptools_dir, "plotCorrelation.tsv")
    log:
        os.path.join(deeptools_dir, "log", "plotCorrelation.log")
    shell:
        "plotCorrelation "
        " -in {input} -c pearson -p heatmap -o {output.corr} "
        " --outFileCorMatrix {output.table} "
        " --colorMap YlOrRd "
        " --removeOutliers"


rule plotPCA:
    input:
        os.path.join(deeptools_dir, "multiBamSummary.results.npz")
    output:
        os.path.join(deeptools_dir, "plotPCA.png")    
    log:
        os.path.join(deeptools_dir, "log", "ploPCA.log")
    shell:
        "plotPCA -in {input} -o {output}"
