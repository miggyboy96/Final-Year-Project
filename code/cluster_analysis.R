## cluster.analysis.R
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
load(file = 'results/output_me.Rda')

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
xlabel <- mapply(paste,colnames(adj_mat),rep("(",length(colnames(adj_mat))),snp_to_gene[colnames(adj_mat)],rep(")",length(colnames(adj_mat))))
grn_effect_pheatmap <- pheatmap(mat = adj_mat, labels_col = xlabel, color = magma(n = 80))

## Reordering snps and genes based on heatmap
connected_snps <- colnames(adj_mat[,grn_effect_pheatmap$tree_col[["order"]]])
regulatory_genes <- rownames(adj_mat[grn_effect_pheatmap$tree_row[['order']],])

## Extract variant and gene clusters based on height of heirarchy
plot(grn_effect_pheatmap$tree_col, main = 'Variant Cluster Dendrogram')
abline(h=30, col="red", lty=2, lwd=2) # height = 30
snp_clusters <- sort(cutree(grn_effect_pheatmap$tree_col, h=30))
snp_cluster_dendrogram <- recordPlot
plot(grn_effect_pheatmap$tree_row, main = 'Regulatory Gene Cluster Dendrogram')
abline(h=40, col="red", lty=2, lwd=2)
gene_clusters <- sort(cutree(grn_effect_pheatmap$tree_row, h=40))

## Write snp and gene lists to .csv files
for (x in 1:max(snp_clusters)){
  snpid <- names(which(snp_clusters==x)) # Creates a vector list of snps in cluster x
  geneid <- snp_to_gene[snpid]
  filename <- paste0("results/clusters/snp_cluster_",x,".csv")
  write.csv(cbind(snpid,geneid),file = filename, row.names = F)
}
for (y in 1:max(gene_clusters)){
  geneid <- names(which(gene_clusters==y))
  filename <- paste0('results/clusters/variant_cluster_',y,'.csv')
  write.csv(geneid, file = filename, row.names = F)
}

## Saving variables
save(list = c('regulatory_genes', 'connected_snps', 'clusters', 'grn_effect_pheatmap', 'cluster_dendrogram'),
     file = 'results/output_cluster.Rda')


