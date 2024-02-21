## 4_heatmap.R
## Miguel Alburo
## 20/02/2023
## cluster analysis of eQTL results

## Packages
library("tidyverse")

## Load data
load(file = 'results/volcano.Rda')
load(file = 'data/processed_data.Rda')
source(file = 'code/functions.R')

# Function to obtain data frame of variant positions and pvalues for a given gene
gene2Manhattan <- function(genename){
  ManhattanData <- res %>%
    inner_join(snpsloc, by = c("SNP" = "snpid")) %>%
    filter(gene == genename) %>%
    select(SNP, chr, pos, FDR) %>%
    rename(P = FDR, CHR = chr, BP = pos)
  return(ManhattanData)
}

# Function to create Manhattan Plot
myManhattan(testy,suggestiveline = FALSE, genomewideline = 10e-6, highlight = intersect(unique(sigres$SNP), unique(testy$SNP)), point.size = 3, label.snps = intersect(unique(sigres$SNP), unique(testy$SNP)))
manplot + geom_label_repel(data=subset(don, is_annotate=="yes"), aes(label=SNP), size=2) +
subse