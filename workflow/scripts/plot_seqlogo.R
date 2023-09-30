library(ggseqlogo)
library(readr)
library(ggplot2)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)


sample <- suppressMessages(read_tsv(args[1], col_names = FALSE))
print(sample)

if (length(sample) > 0){
    # print("YESSS")
    png(args[2], width=800, height=200)
    print(ggseqlogo(sample$X1, seq_type="rna"))
    dev.off()
} else {
    # print("NOOO")
    png(args[2])
    print(ggplot())
    dev.off()
}
