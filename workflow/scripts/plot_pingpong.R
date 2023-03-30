library(dplyr)
library(readr)
library(ggplot2)
library(ggridges)
library(stringr)

read_mutate_tsv <- function(path){
    tab <- read_tsv(path)
    i <- path
    filename <- str_split(i, "/")[[1]]
    basename <- str_split(filename[length(filename)], "\\.")[[1]][1]
    tab <- tab %>% mutate(sample = basename)
    return(tab)
}


# pingpong_files <- c(
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg1.piRNAs.pingpong.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg2.piRNAs.pingpong.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym1.piRNAs.pingpong.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym2.piRNAs.pingpong.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs1.piRNAs.pingpong.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs2.piRNAs.pingpong.tsv"
# )

pingpong_files <- snakemake@input
ridges <- snakemake@output[[1]]
lines <- snakemake@output[[2]]

pingpongs <- lapply(pingpong_files, read_mutate_tsv)
pingpongs_z <- lapply(pingpongs, function(x){
    zs <- x$count
    total <- sum(x$count)
    x$zscore <- zs
    pct <- sapply(x$count, function(y){return(y * 100 / total)})
    x$pct <- pct
    return(x)
})

ppz_bind <- bind_rows(pingpongs_z)
ppz_bind$sample <- factor(ppz_bind$sample, levels = c("nym2", "nym1",  "emb2", "emb1", "egg2", "egg1", "pvs2", "pvs1"))

png(ridges, width = 800)
ggplot(ppz_bind, aes(overlap_length, sample, height = pct, group = sample, fill = sample)) + 
  geom_density_ridges(stat="identity", scale = 4, alpha = 0.5, show.legend = FALSE) +
  scale_x_continuous(breaks = c(1,5,10,15,20), limits = c(1,20)) +
  theme_ridges(font_size = 18) +
  ylab("") + xlab("Overlap Length (nt)")
dev.off()

ppz_bind$sample <- factor(ppz_bind$sample, levels = c("pvs1", "pvs2",  "egg1", "egg2", "emb1", "emb2", "nym1", "nym2"))
png(lines, width = 800)
ggplot(ppz_bind, aes(x = overlap_length, y = pct, color=sample)) + 
    geom_line(linewidth=1) +
    theme_bw(base_size = 20) +
    theme(legend.title = element_blank()) +
    scale_x_continuous(breaks = c(1,10,20), limits = c(1,20)) +
    ylab("Pairs (%)") + xlab("Overlap Length (nt)")
dev.off()
