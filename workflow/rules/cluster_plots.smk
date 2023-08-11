
cluster_plots = os.path.join(cluster_dir, "plots")

rule plot_ideograms:
    input:
        clusters = os.path.join(cluster_dir, "clusters.bb"),
        pvs = os.path.join(cluster_dir, "pvs.clusters.nomir.bed"),
        egg = os.path.join(cluster_dir, "egg.clusters.nomir.bed"),
        emb = os.path.join(cluster_dir, "emb.clusters.nomir.bed"),
        nym = os.path.join(cluster_dir, "nym.clusters.nomir.bed")
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

rule plot_clusters_heatmap:
    input:
        outcount = os.path.join(cluster_dir, "clusters.counts.tsv"),
        final = os.path.join(cluster_dir, "clusters.final.tsv")
    output:
        os.path.join(cluster_plots, "clusters.heatmap.zscore.png")
    script:
        "../scripts/plot_clusters_heatmap.R"


