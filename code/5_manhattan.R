## expression_heatmap.R
## Miguel Alburo
## 20/02/2024
## cluster analysis of eQTL results

## Packages
library("tidyverse")
library(RColorBrewer)

## Load data
load(file = 'results/volcano.Rda')
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')
source(file = '1_functions.R')

# Function to obtain data frame of variant positions and pvalues for a given gene
gene2Manhattan <- function(genename){
  ManhattanData <- output_eqtl %>%
    inner_join(snpsloc, by = c("SNP" = "snpid")) %>%
    filter(gene == genename) %>%
    select(SNP, chr, pos, FDR) %>%
    rename(P = FDR, CHR = chr, BP = pos)
  return(ManhattanData)
}

manplots <- sapply(unique(sigres$gene),
                   function(g) {
                     longdata <- gene2Manhattan(g)
                     myManhattan(longdata,graph.title = g,col = brewer.pal(8, "Dark2"),suggestiveline = FALSE, genomewideline = 10e-6, highlight = intersect(unique(sigres$SNP), unique(longdata$SNP)), point.size = 3, label.snps = intersect(unique(sigres$SNP), unique(longdata$SNP)))
                   }, simplify = FALSE)

for (plt in seq_along(manplots)){
  print(manplots[plt])
}