## 4_expression_heatmap.R
## Miguel Alburo
## 20/02/2023
## cluster analysis of eQTL results

## Packages
library("tidyverse")
library(RColorBrewer)

## Load data
load(file = 'results/volcano.Rda')
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')
source(file = 'code/functions.R')

# Function to obtain data frame of variant positions and pvalues for a given gene
gene2Manhattan <- function(genename){
  ManhattanData <- output_eqtl %>%
    inner_join(snpsloc, by = c("SNP" = "snpid")) %>%
    filter(gene == genename) %>%
    select(SNP, chr, pos, FDR) %>%
    rename(P = FDR, CHR = chr, BP = pos)
  return(ManhattanData)
}

testy <- gene2Manhattan("AT4G32800")
# Function to create Manhattan Plot
myManhattan(testy,col = brewer.pal(8, "Dark2"),suggestiveline = FALSE, genomewideline = 10e-6, highlight = intersect(unique(sigres$SNP), unique(testy$SNP)), point.size = 3, label.snps = intersect(unique(sigres$SNP), unique(testy$SNP)))
#

myManhattan(testy, highlight = 1e-6, point.size = 3)