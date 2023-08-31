library(ggseqlogo)
library(readr)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)


sample <- suppressMessages(read_tsv(args[1], col_names = FALSE))

png(args[2], width=800, height=200)
ggseqlogo(sample$X1, seq_type="rna")
dev.off()
