## 4_expression_heatmap.R
## Miguel Alburo
## 20/02/2023
## cluster analysis of eQTL results

## Load packages
library(readxl) # Reading excel spreadsheets .xlsx
library(dplyr)  # Data manipulation tools.
library(pheatmap) # Heatmap vizualisation.
library(viridis)  # Heatmap visualization, colors.
library(RColorBrewer) # Heatmap colours

## Loading data
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')
load(file = 'results/volcano.Rda')

# List of variants and genes
snplist <- unique(sigres$SNP)
genelist <- unique(sigres$gene)

## Constructing an adjacency matrix
adj_mat <- matrix(
  data = 0,
  dimnames = list(genelist, snplist),
  nrow = length(genelist), ncol = length(snplist)
)
for (i in seq(nrow(sigres))){
  snp <- sigres$SNP[i]
  gene <- sigres$gene[i]
  pi <- sigres$log2FC[i]
  adj_mat[gene,snp] <- pi
}

## Color palette
cool <- rainbow(50, start=rgb2hsv(col2rgb('cyan'))[1], end=rgb2hsv(col2rgb('blue'))[1])
warm <- rainbow(50, start=rgb2hsv(col2rgb('red'))[1], end=rgb2hsv(col2rgb('yellow'))[1])
cols <- c(rev(cool), "snow2", "snow2", rev(warm))
mypalette <- colorRampPalette(cols)(255)

## Printing the heatmap
xlabel <- mapply(paste,colnames(adj_mat),rep("(",length(colnames(adj_mat))),snp_to_gene[colnames(adj_mat)],rep(")",length(colnames(adj_mat))))
grn_effect_pheatmap <- pheatmap(mat = adj_mat, labels_col = xlabel, color = mypalette)

## Reordering snps and genes based on heatmap
sigvariants <- colnames(adj_mat[,grn_effect_pheatmap$tree_col[["order"]]])
regulatorygenes <- rownames(adj_mat[grn_effect_pheatmap$tree_row[['order']],])

## Extract variant and gene clusters based on height of heirarchy
plot(grn_effect_pheatmap$tree_col, main = 'Variant Cluster Dendrogram')
abline(h=30, col="red", lty=2, lwd=2) # height = 30
snp_clusters <- sort(cutree(grn_effect_pheatmap$tree_col, h=30))
snp_cluster_dendrogram <- recordPlot()
plot(grn_effect_pheatmap$tree_row, main = 'Regulatory Gene Cluster Dendrogram')
abline(h=40, col="red", lty=2, lwd=2)
gene_clusters <- sort(cutree(grn_effect_pheatmap$tree_row, h=40))
gene_cluster_dendrogram <- recordPlot()

## Write snp and gene lists to .csv files
for (x in 1:max(snp_clusters)){
  snpid <- names(which(snp_clusters==x)) # Creates a vector list of snps in cluster x
  geneid <- snp_to_gene[snpid]
  filename <- paste0("results/clusters/snp_cluster_",x,".csv")
  write.csv(cbind(snpid,geneid),file = filename, row.names = F)
}
for (y in 1:max(gene_clusters)){
  geneid <- names(which(gene_clusters==y))
  filename <- paste0('results/clusters/gene_cluster_',y,'.csv')
  write.csv(geneid, file = filename, row.names = F)
}

## Saving variables
save(list = c('gene_clusters', 'snp_clusters', 'grn_effect_pheatmap', 'gene_cluster_dendrogram', 'snp_cluster_dendrogram'),
     file = 'results/clusters/output_cluster.Rda')


