library(readr)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggridges)
library(forcats)

read_mutate_tsv <- function(path){
    tab <- read_tsv(path)
    i <- path
    filename <- str_split(i, "/")[[1]]
    basename <- str_split(filename[length(filename)], "\\.")[[1]][1]
    tab <- tab %>% mutate(sample = basename)
    return(tab)
}

# print(snakemake@input)

# files <- c(
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/230120_01_nymph_rhodnius_rep1_R1.sorted.lengths.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/230120_10_nymph_rhodnius_rep2_R1.sorted.lengths.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/230120_04_embryo_rhodnius_rep2_R1.sorted.lengths.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/230120_09_embryo_rhodnius_rep1_R1.sorted.lengths.tsv"
# )
# tables <- lapply(files, read_mutate_tsv)

tables <- lapply(snakemake@input, read_mutate_tsv)

merged_table <- bind_rows(tables)

png(snakemake@output[[1]], height = 400, width = 1000)
ggplot(merged_table, aes(length, sample, height = total, group = sample, fill = sample)) + 
#   geom_ridgeline(scale=1) 
  geom_density_ridges(stat="identity", scale = 4, alpha = 0.5) +
  scale_x_continuous(breaks = c(18:35), limits = c(18,35)) +
  theme_ridges() +
  ylab("")
#   xlim(18, 35)
dev.off()