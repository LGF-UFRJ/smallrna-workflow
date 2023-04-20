
cluster_plots = os.path.join(cluster_dir, "plots")

rule plot_ideograms:
    input:
        pvs = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/pvs.clusters.bed",
        egg = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/egg.clusters.bed",
        emb = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/emb.clusters.bed",
        nym = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/nym.clusters.bed"
    output:
        os.path.join(cluster_plots, "clusters.karyoplot.png")
    params:
        config["genome_index"]
    script:
        "../scripts/karyoplot.R"

rule plot_cluster_size:
    input:
        os.path.join(cluster_dir, "clusters.final.tsv")
    output:
        os.path.join(cluster_plots, "clusters.sizes.png")
    script:
        "../scripts/cluster_size_chr.R"

rule plot_cluster_upset:
    input:
        os.path.join(cluster_dir, "clusters.overlap.tsv")
    output:
        os.path.join(cluster_plots, "clusters.overlap.upset.png")
    script:
        "../scripts/cluster_upset.R"
