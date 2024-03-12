## 3_volcano_plot.R
## Miguel Alburo
## 19/02/2024
## Volcano plot from MatrixEQTL results

## Packages
library(tidyverse)
library(EnhancedVolcano)

## Loading data
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')
source(file = 'code/functions.R')

## Gene selection
res <- output_eqtl %>%
  filter(gene %in% regulatorygenes) %>%
  select(SNP, gene, FDR) %>%  # Adjust these column names as per your data frame
  rowwise() %>%
  mutate(log2FC = calclog2FC(SNP,gene)) %>%
  mutate(pi = -1 * log10(FDR) * calclog2FC(SNP,gene)) %>%
  ungroup()

# Volcano Plot
EnhancedVolcano(res,
  lab = res$SNP,
  x = 'log2FC',
  y = 'FDR',
  title = 'MatrixEQTL results',
  pCutoff = 0.05,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 6.0)

# Subset res based on volcano thresholds
sigres <- res %>%
  filter((log2FC>1 | log2FC < -1) & FDR < 0.05)
sigvariants <- unique(sigres$SNP)
sig_genes <- unique(sigres$gene)

# Write genelist to csv
snpgenes <- snp_to_gene[sigvariants]
write.csv(snpgenes, file = "results/sigvariants.csv", row.names = T)

# Print results
print(paste(length(unique(sigres$SNP)),"significant variants")) # 58 variants
print(paste(length(unique(sigres$gene)),"differentially expressed regulators")) # 19 variants

# Save variables
save(list = c('res', 'sigres', 'sigvariants', 'sig_genes', 'snpgenes'), file = "results/volcano.Rda")
rm(list=ls())