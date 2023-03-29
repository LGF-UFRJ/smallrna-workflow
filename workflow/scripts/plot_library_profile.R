library(readr)
library(rtracklayer)
library(dplyr)
library(ggplot2)
library(stringr)
library(tidyr)

get_basenames <- function(map_count_files){
    split_folders <- sapply(map_count_files, function(x){return(str_split(x, "/"))})
    remove_ext <- sapply(split_folders, function(x){return(str_split(x[length(x)], "\\."))})
    basename_vec <- sapply(remove_ext, function(x){return(x[1])})
    names(basename_vec) <- NULL
    return(basename_vec)
}

# Inputs
# count_file <- "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/library_profile/counts.tsv"
# gff_file <- "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/library_profile/VectorBase-61_RprolixusCDC.gff"
# map_count_files <- c(
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs1.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/pvs2.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg1.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/egg2.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb1.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/emb2.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym1.sorted.counts.tsv",
#     "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/results/mapping/nym2.sorted.counts.tsv"
# )

map_count_files <- snakemake@input
gff_file <- snakemake@params[[1]]
count_file <- snakemake@params[[2]]
output <- snakemake@output[[1]]

# Importing datasets
counts <- read_tsv(count_file)
gff <- import(gff_file)

# Getting unassigned reads
## Get total mapped reads for each sample
names(map_count_files) <- get_basenames(map_count_files)
map_count <- lapply(map_count_files, function(x){
    return(read_tsv(x, col_names = FALSE))
})

## Get total assigned reads for each sample
samples <- colnames(counts)[c(7:length(colnames(counts)))]
names(samples) <- samples
total_assigned <- lapply(samples, function(x){
    return(sum(counts[x]))
})

## Get unassigned (total mapped - assigned)
unassigned <- list()
for (sample in names(total_assigned)){
    total <- map_count[[sample]][1,2] %>% as.numeric
    assigned <- total_assigned[[sample]]
    unassigned[[sample]] <- total - assigned
}

# Characterizing library
gff_df <- as.data.frame(gff)
names(gff_df)[10] <- "Geneid"

counts_gff <- left_join(counts, gff_df, by = "Geneid") %>%
    mutate(
        biotype = case_when(
            grepl("rRNA", Geneid) | grepl("rRNA", ebi_biotype) ~ "rRNA", 
            # grepl("snRNA", Geneid) | grepl("snRNA", ebi_biotype) ~ "snRNA", 
            grepl("tRNA", Geneid) | grepl("tRNA", ebi_biotype) ~ "tRNA",
            grepl("miRNA", ebi_biotype) ~ "miRNA",
            is.na(ebi_biotype) & !grepl("rRNA|snRNA|tRNA", Geneid) ~ "repeat",
            grepl("protein_coding", ebi_biotype) ~ "gene",
            TRUE ~ "other"
        )
    )

counts_gff_sel <- counts_gff[c(1,7:14,25,32)]
counts_sum <- counts_gff_sel %>%
    group_by(biotype) %>%
    summarise(pvs1 = sum(pvs1), pvs2 = sum(pvs2),
              egg1 = sum(egg1), egg2 = sum(egg2), 
              emb1 = sum(emb1), emb2 = sum(emb2), 
              nym1 = sum(nym1), nym2 = sum(nym2))


unassigned_row <- c("biotype" = "unassigned")
for (sample in samples){
    unassigned_row <- c(unassigned_row, unassigned[sample])
}
unassigned_row

# Calculating percentage
counts_pct <- bind_rows(counts_sum, unassigned_row)
for (sample in samples){
    sample_count_sum <- counts_pct[sample] %>% pull
    total_mapped <- map_count[[sample]][1,2] %>% as.numeric
    sample_pct <- sapply(sample_count_sum, function(x){return(x * 100 / total_mapped)})
    counts_pct[sample] <- sample_pct
}
counts_pct

counts_pct_long <- pivot_longer(counts_pct, c("pvs1", "pvs2", "egg1", "egg2", "emb1", "emb2", "nym1", "nym2"))

counts_pct_long$name <- factor(counts_pct_long$name, levels = c("nym2", "nym1",  "emb2", "emb1", "egg2", "egg1", "pvs2", "pvs1"))

counts_pct_long$name <- factor(counts_pct_long$name, levels = c("nym2", "nym1",  "emb2", "emb1", "egg2", "egg1", "pvs2", "pvs1"))

colors <- c("#006BA6", "#F9DC5C" ,"#A1C084" ,"#F26430" ,"#A846A0" ,"grey", "grey25")

png(output)
ggplot(counts_pct_long, mapping = aes(y = name, x = value, fill = biotype)) +
    geom_col() +
    # scale_fill_brewer(palette = "Dark2") +
    scale_fill_manual(values = c("#006BA6", "#F9DC5C", "grey","#F26430" ,"#A846A0" ,"#A1C084", "grey25")) +
    theme_bw(base_size = 18) + 
    ylab("") + xlab("Percentage (%)") + labs(fill = "")
dev.off()

