## 3_volcano_plot.R
## Miguel Alburo
## 19/02/2024
## Volcano plot from MatrixEQTL results

## Packages
library(ggrepel)
library(EnhancedVolcano)
library(readxl)
library(dplyr)

## Loading data
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')
source(file = 'code/functions.R')

## Gene selection
# Importing GRN data
gene_regulatory_network <- read_excel(path = 'data/raw/media-10.xlsx') # Downloaded from bioRxiv's Supplementary Materials
regulatorygenes <- unique(gene_regulatory_network$`regulatory gene`)
res <- output_eqtl[which(output_eqtl$gene %in% regulatorygenes),]

## Function to calculate effect size as log2FC
calclog2FC <- function(snp, gene){
  # Extract expression values for the gene
  gene_expr <- expr[expr$geneid == gene, ]
  # Extract genotype values for the SNP
  variant_genotype <- gt[gt$snpid == snp, ]
  # Ensure matching samples between expression and genotype data
  common_samples <- intersect(names(gene_expr), names(variant_genotype))
  gene_expr <- gene_expr[common_samples]
  variant_genotype <- variant_genotype[common_samples]
   # Group expression values by genotype (0 vs 1 or 2)
  expr_genotype_before <- gene_expr[variant_genotype == min(variant_genotype)]
  expr_genotype_after <- gene_expr[variant_genotype == max(variant_genotype)] # either 1 or 2
  # Calculate average expression for each genotype group
  avg_expr_genotype_before <- mean(expr_genotype_before, na.rm = TRUE)
  avg_expr_genotype_after <- mean(expr_genotype_after, na.rm = TRUE)
  # Calculate log2 fold change
  log2fc <- log2(avg_expr_genotype_after / avg_expr_genotype_before)
  return(log2fc)
}

## Data frame for plot ~ 30 seconds
res <- res %>%
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
  pCutoff = 10e-6,
  FCcutoff = 1,
  pointSize = 3.0,
  labSize = 6.0)

# Subset res based on volcano thresholds
sigres <- res %>%
  filter((log2FC>1 | log2FC < -1) & FDR < 10e-6)
length(unique(sigres$SNP)) # 76 variants
length(unique(sigres$gene)) # 20 genes
# Save variables
save(list = c('res', 'sigres', 'regulatorygenes'), file = "results/volcano.Rda")
rm(list=ls())