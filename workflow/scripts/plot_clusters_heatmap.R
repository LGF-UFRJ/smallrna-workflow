library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)
library(pheatmap)
library(tibble)
library(RColorBrewer)

counts <- read_tsv("/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/clusters/clusters.counts.tsv")
n <- names(counts)
n <- gsub("#", "", n)
n <- gsub("\'", "", n)
n <- gsub(".sorted.bam", "", n)
names(counts) <- n

final <- read_tsv("/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/clusters/clusters.final.tsv")

nmapped <- list(
    "pvs1" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs1.sorted.counts.tsv",
    "pvs2" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs2.sorted.counts.tsv",
    "egg1" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg1.sorted.counts.tsv",
    "egg2" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg2.sorted.counts.tsv",
    "emb1" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb1.sorted.counts.tsv",
    "emb2" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb2.sorted.counts.tsv",
    "nym1" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym1.sorted.counts.tsv",
    "nym2" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym2.sorted.counts.tsv"
)

nmapped_total <- lapply(nmapped, function(x){
    tab <- read_tsv(x, col_names = FALSE)
    tab <- pivot_wider(tab, names_from = X1, values_from = X2)
    return(tab$total)
})

samples <- names(counts)[4:length(names(counts))]

for (i in samples){
    total <- nmapped_total[[i]]
    sample_counts <- counts[[i]]
    norm_counts <- (sample_counts * 1000000) / total
    counts[[i]] <- norm_counts
}

final_counts <- left_join(final, counts, by = c("chr", "start", "end"))

mat <- column_to_rownames(final_counts, var = "name") %>% select(all_of(samples))

palette <- colorRampPalette(brewer.pal(11,'YlOrRd'))
png("plot.png")
pheatmap(mat[c(1:30),], cluster_cols = FALSE, scale="row", angle_col = 45)
dev.off()
