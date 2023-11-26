## grn_effect.R
## Miguel Alburo
## 26/11/2023
## Analysis of significant eQTLs with high connection to regulatory genes in the GRN.

## Load packages
  library(readxl) # Reading excel spreadsheets .xlsx
  library(dplyr)  # Data manipulation tools.
  library(pheatmap) # Heatmap vizualisation.
  library(viridis)  # Heatmap visualization, colors.

## Loading data
  load(file = 'data/processed_data.Rda')
  load(file = 'data/output_me.Rda')

## Importing GRN data
  gene_regulatory_network <- read_excel(path = 'data/raw/media-10.xlsx') # Downloaded from bioRxiv's Supplementary Materials

## Regulatory Genes
  regulatory_genes <- unique(gene_regulatory_network$`regulatory gene`)

## Variants which are highly connected to regulatory genes.
  subset_output_eqtl <- output_eqtl[which(output_eqtl$gene %in% regulatory_genes),] # Subset df for regulatory genes
  snp_table<- table(subset_output_eqtl$SNP) # Tables number of regulatory genes each snp is associated with
  snp_table <- snp_table[order(-snp_table)] # Descending order
  connected_snps <- unlist(dimnames(snp_table[1:100])) # 100 most connected snps
  subset_output_eqtl <- subset_output_eqtl[which(subset_output_eqtl$SNP %in% connected_snps),] # Subset for connected snps

## Constructing an adjacency matrix
  adj_mat <- matrix(
    data = 0,
    dimnames = list(regulatory_genes, connected_snps),
    nrow = length(regulatory_genes), ncol = length(connected_snps)
  )
  for (i in seq(nrow(subset_output_eqtl))){
    snp <- subset_output_eqtl$SNP[i]
    gene <- subset_output_eqtl$gene[i]
    logp <- -log(subset_output_eqtl$p.value[i], base = 10)
    adj_mat[gene,snp] <- logp
  }

## Printing the heatmap
  grn_effect_pheatmap <- pheatmap(mat = adj_mat, color = magma(n = 80))