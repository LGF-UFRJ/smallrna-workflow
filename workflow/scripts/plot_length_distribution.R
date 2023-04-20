library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggridges)
library(forcats)
library(scales)

read_mutate_tsv <- function(path){
    tab <- read_tsv(path)
    i <- path
    filename <- str_split(i, "/")[[1]]
    basename <- str_split(filename[length(filename)], "\\.")[[1]][1]
    tab <- tab %>% mutate(sample = basename)
    return(tab)
}

# files <- list(
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs1.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs2.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg1.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg2.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb1.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb2.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym1.sorted.lengths.tsv",
#   "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym2.sorted.lengths.tsv"
# )

# tables <- lapply(files, read_mutate_tsv)

cols <- c("#81488f", "#cf05c5", "#414cb6", "#3990d8", "#0c9673", "#38792b", "#b9b603", "#ca4a4a")

tables <- lapply(snakemake@input, read_mutate_tsv)

merged_table <- bind_rows(tables)
merged_table$sample <- factor(merged_table$sample, levels = c("nym2", "nym1",  "emb2", "emb1", "egg2", "egg1", "pvs2", "pvs1"))

png(snakemake@output[[1]], height = 400, width = 1000)
ggplot(merged_table, aes(length, sample, height = total, group = sample, fill = sample)) + 
  geom_density_ridges(stat="identity", scale = 4, alpha = 0.5, show.legend = FALSE) +
  scale_x_continuous(breaks = c(18:35), limits = c(18,30)) +
  theme_ridges(font_size = 20) +
  ylab("") + xlab("Length (nt)")
dev.off()



merged_table$sample <- factor(merged_table$sample, levels = c("pvs1", "pvs2",  "egg1", "egg2", "emb1", "emb2", "nym1", "nym2"))
png(snakemake@output[[2]], width = 800)
ggplot(merged_table, aes(x = length, y = total, color=sample)) + 
    scale_shape_manual(values = c(1:8)) +
    scale_color_manual(values = cols) +
    scale_y_continuous(labels = comma) +
    geom_line(linewidth=1) +
    scale_x_continuous(breaks = c(18:35), limits = c(18,30)) +
    theme_bw(base_size = 20) +
    ylab("RPM") + xlab("Length (nt)")
dev.off()
