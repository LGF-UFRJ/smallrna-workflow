library(dplyr)
library(readr)
library(UpSetR)


overlap_file <- snakemake@input[[1]]
outfile <- snakemake@output[[1]]

overlap <- read_tsv(overlap_file)


overlap_fmt <- as.data.frame(overlap)

samples <- c("pvs1", "pvs2", "egg1", "egg2", "emb1", "emb2", "nym1", "nym2")
for (i in samples){
    overlap_fmt[i] <- ifelse(overlap_fmt[i] == 0, 0, 1)
}

png(outfile, width = 800)
upset(overlap_fmt, 
    sets = rev(samples), 
    order.by = "freq",
    # keep.order = TRUE,
    sets.x.label = "Number of piRNA Clusters",
    mainbar.y.label = "piRNA Cluster Intersections",
    text.scale = c(1.5, 1.5, 1.5, 1.3, 1.5, 1.3))
dev.off()
