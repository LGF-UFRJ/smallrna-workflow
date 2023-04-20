library(dplyr)
library(readr)
library(ggplot2)

final_file <- snakemake@input[[1]]

outfile <- snakemake@output[[1]]

final <- read_tsv(final_file)

chrs <- c()
for (i in seq(11)){
    c <- paste0("HiC_scaffold_", i)
    chrs <- c(chrs, c)
}

final_fmt <- final %>%
    filter(chr %in% chrs)
final_fmt$chr <- factor(final_fmt$chr, levels = chrs)


png(outfile)
ggplot(final_fmt, mapping = aes(x = chr, y = length, fill = chr)) +
    geom_boxplot(show.legend = FALSE) +
    scale_y_log10(
        breaks = scales::trans_breaks("log10", function(x) 10^x),
        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_bw(base_size = 16) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size=16, color = "black")) +
    xlab("") + ylab("Length (bp)")
dev.off()


