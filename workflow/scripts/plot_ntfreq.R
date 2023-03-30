library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(ggridges)
library(tidyr)

read_mutate_tsv <- function(path){
    tab <- read_tsv(path)
    i <- path
    filename <- str_split(i, "/")[[1]]
    basename <- str_split(filename[length(filename)], "\\.")[[1]][1]
    tab <- tab %>% mutate(sample = basename)
    return(tab)
}


# files <- c(
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs1.piRNAs.sense.ntfreq.tsv",
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs2.piRNAs.sense.ntfreq.tsv",
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg1.piRNAs.sense.ntfreq.tsv",
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg2.piRNAs.sense.ntfreq.tsv",
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym1.piRNAs.sense.ntfreq.tsv",
# #     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym2.piRNAs.sense.ntfreq.tsv"
# # )

#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs1.piRNAs.antisense.ntfreq.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/pvs2.piRNAs.antisense.ntfreq.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg1.piRNAs.antisense.ntfreq.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/egg2.piRNAs.antisense.ntfreq.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym1.piRNAs.antisense.ntfreq.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/piRNAs/nym2.piRNAs.antisense.ntfreq.tsv"
# )

files <- snakemake@input
outfile <- snakemake@output[[1]]

ntfreq <- lapply(files, read_mutate_tsv)

ntfreq_bind <- bind_rows(ntfreq)

ntfreq_long <- pivot_longer(ntfreq_bind, c("T", "A", "C", "G", "N"))

ntfreq_long$name <- factor(ntfreq_long$name, levels = c("T", "A", "C", "G", "N"))
ntfreq_long$sample <- factor(ntfreq_long$sample, levels = c("pvs1", "pvs2",  "egg1", "egg2", "emb1", "emb2", "nym1", "nym2"))
png(outfile, width = 1400, height = 400)
ggplot(ntfreq_long %>% filter(name != "N") %>% filter(sample != "egg1"), aes(x = position, y = value, color=sample)) + 
    geom_line(linewidth=1) +
    theme_bw(base_size = 30) +
    theme(legend.title = element_blank()) +
    facet_wrap(vars(name), nrow = 1) +
    scale_x_continuous(breaks = c(1,10,20,30), limits = c(1,30)) +
    ylab("Frequency (%)") + xlab("Position")
dev.off()

