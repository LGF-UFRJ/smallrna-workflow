
rule get_miRNAs:
    input:
        bam = lambda wildcards: os.path.join(featureCount_dir, wildcards.sample + ".sorted.nh.bam.featureCounts.bam"),
        names = os.path.join(featureCount_dir, "annotations", "names.tsv")
    output:
        os.path.join(mirna_dir, "{sample}.miRNAs.bam")
    log:
        os.path.join(map_out_vb_dir, "log", "{sample}.index_uniquely.log")
    #conda:
        #"../envs/pysam0.yaml"
    script:
        "../scripts/get_miRNAs.py"
