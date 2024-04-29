## 3_expression_heatmap.R
## Miguel Alburo
## 20/02/2024
## cluster analysis of eQTL results

## Load packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer) # Heatmap colours

## Loading data
load(file = 'data/processed_data.Rda')
load(file = 'results/output_me.Rda')

# List of variants and genes
snplist <- sigvariants
genelist <- sig_genes

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

# Column adjustments
#filtered_snps <- names(sort(table(sigres$SNP), decreasing = T)[1:100])
#adj_mat <- adj_mat[,filtered_snps]
colnames(adj_mat) <- mapply(paste0,colnames(adj_mat),rep(" (",length(colnames(adj_mat))),snp_to_gene[colnames(adj_mat)],rep(")",length(colnames(adj_mat))))

## Color palette
col_fun <- colorRamp2(c(min(adj_mat, na.rm = T), -0.5, 0, 0.5, max(adj_mat, na.rm = T)),
                      c("blue", "cyan", "snow3", "yellow","red")
)

## Printing the heatmap
expression_heatmap <- Heatmap(mat = adj_mat,
        name = "log2FC",
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 0.1),
        row_title = "Regulatory genes",
        row_title_side = "left",
        row_names_side = "left",
        column_title = "Variants",
        column_title_side = "top",
        show_column_names = F
)
                  # labels_col = NULL,
                  # show_colnames = F,
                  # col = mypalette,
                  # breaks = c(bk1,0,bk2))

gene_table <- as.data.frame(sort(table(sigres$gene), decreasing = T))
# FILL
colors <- c(rep("> 40", 11), "> 20", "> 20", rep("0 - 20", 7)) # Since 13 significant regulatory genes


p.gene_table <- ggplot(gene_table, aes(x=Var1, y=Freq, fill = colors))+
  geom_bar(stat="identity") +
  labs(x="Bolting regulator", y="Number of variants") +
  scale_fill_manual(values=c("orange", "red", "darkgrey")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## Reordering snps and genes based on heatmap
# sigvariants <- colnames(adj_mat[,expression_heatmap$tree_col[["order"]]])
# sig_genes <- rownames(adj_mat[expression_heatmap$tree_row[['order']],])

## Extract variant and gene clusters based on height of heirarchy
# plot(expression_heatmap$tree_col, main = 'Variant Cluster Dendrogram')
# abline(h=30, col="red", lty=2, lwd=2) # height = 30
# snp_clusters <- sort(cutree(expression_heatmap$tree_col, h=30))
# snp_cluster_dendrogram <- recordPlot()
# plot(expression_heatmap$tree_row, main = 'Regulatory Gene Cluster Dendrogram')
# abline(h=40, col="red", lty=2, lwd=2)
# gene_clusters <- sort(cutree(expression_heatmap$tree_row, h=40))
# gene_cluster_dendrogram <- recordPlot()

## Write snp and gene lists to .csv files
# for (x in 1:max(snp_clusters)){
#   snpid <- names(which(snp_clusters==x)) # Creates a vector list of snps in cluster x
#   geneid <- snp_to_gene[snpid]
#   filename <- paste0("results/clusters/snp_cluster_",x,".csv")
#   write.csv(cbind(snpid,geneid),file = filename, row.names = F)
# }
# for (y in 1:max(gene_clusters)){
#   geneid <- names(which(gene_clusters==y))
#   filename <- paste0('results/clusters/gene_cluster_',y,'.csv')
#   write.csv(geneid, file = filename, row.names = F)
# }

# Print results


## Saving variables
save(list = c('gene_clusters', 'snp_clusters', 'expression_heatmap', 'gene_cluster_dendrogram', 'snp_cluster_dendrogram'),
     file = 'results/clusters/output_cluster.Rda')

ggsave(p.gene_table, file ="results/plots/gene_table.pdf", width=90, height=90, units = "mm", dpi=300)

pdf(file = "results/plots/expression_heatmap.pdf", width = 190 * 0.0393701, height = 110 * 0.0393701)
draw(expression_heatmap)
dev.off()

rm(list=ls())

print()