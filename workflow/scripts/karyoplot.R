library(karyoploteR)
library(rtracklayer)
library(readr)
library(dplyr)

fai <- read_tsv(snakemake@params[[1]], col_names = FALSE)
fai <- fai %>%  mutate(s=(0*X2)+1) %>% select(X1,s,X2)

custom_genome <- toGRanges(fai %>% as.data.frame)
custom_genome[c(1:11),]

# cluster_files <- c(
#     "pvs" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/pvs.clusters.bed",
#     "egg" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/egg.clusters.bed",
#     "emb" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/emb.clusters.bed",
#     "nym" = "/data/iovino/group2/brito/analyses/smallrna_rhodnius_2023/sandbox/pirnas/nym.clusters.bed"
# )

cluster_files <- snakemake@input

clusters <- lapply(cluster_files, import)

plot.params <- getDefaultPlotParams(plot.type=6)
plot.params$ideogramheight <- 100
plot.params$data1outmargin <- 100


u <- 0.25

bg <- 0 + (4 * u)
p <- 0 + (3 * u) 
g <- 0 + (2 * u)
m <- 0 + (1 * u)
n <- 0 + (0 * u)

png(snakemake@output[[1]], width=1000, height = 800)
kp <- plotKaryotype(genome = custom_genome, chromosomes = fai$X1[c(1:11)], plot.type = 6, plot.params=plot.params)
kpDataBackground(kp, color = "grey90")
kpAbline(kp, h=p+u)
kpPlotRegions(kp, data=clusters[["pvs"]], r0=p, r1=p+u, col = "purple")
kpAbline(kp, h=p)
kpPlotRegions(kp, data=clusters[["egg"]], r0=g, r1=g+u, col = "blue")
kpAbline(kp, h=g)
kpPlotRegions(kp, data=clusters[["emb"]], r0=m, r1=m+u, col = "darkgreen")
kpAbline(kp, h=m)
kpPlotRegions(kp, data=clusters[["nym"]], r0=n, r1=n+u, col = "red")
kpAbline(kp, h=n)
kpAddBaseNumbers(kp, tick.dist = 10000000, tick.len = 10, tick.col="red", cex=1,
                 minor.tick.dist = 1000000, minor.tick.len = 5, minor.tick.col = "black")
dev.off()