## expression_heatmap.R
## Miguel Alburo
## 20/02/2024
## cluster analysis of eQTL results

## Load packages
library(tidyverse)
source(file = "code/1_functions.R")
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer) # Heatmap colours
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

## Color palette
col_fun <- colorRamp2(c(min(adj_mat, na.rm = T), -0.5, 0, 0.5, max(adj_mat, na.rm = T)),
                      c("blue", "cyan", "snow3", "yellow","red")
)

# Kmean clustering
gene_clusters <- kmeans(adj_mat, centers = 4, iter.max = 30, nstart = 2)$cluster
gene_clusters[gene_clusters==gene_clusters["AT3G11440"]] <- gene_clusters["AT2G38470"]
variant_clusters <- kmeans(t(adj_mat), centers = 3, iter.max = 30, nstart = 2)$cluster

# gene_clusters.this <- gene_clusters
# variant_clusters.this <- variant_clusters


## Printing the heatmap
expression_heatmap <- Heatmap(mat = adj_mat,
        name = "log2FC",
        col = col_fun,
        rect_gp = gpar(col = "white", lwd = 0.1),
        row_title_side = "right",
                              row_title_rot = 0,
        row_names_side = "left",
        column_title_side = "top",
        show_column_names = F,
                              column_split = variant_clusters.this,
                              column_title = c("faster bolt", "", "delayed bolt"),
                              row_split =  gene_clusters.this,
                              row_title = c("decreasing TFs", "", "increasing TFs")
)

expression_heatmap <- draw(expression_heatmap, row_title_side = "left", row_title = "Regulatory genes", column_title = "Variants",)


## Connectivity of SNPs and genes
gene_table <- as.data.frame(sort(table(sigres$gene), decreasing = T))
# FILL
group <- c(rep("> 40", 11), "> 20", "> 20", rep("0 - 20", 7)) # Since 13 significant regulatory genes
p.gene_table <- ggplot(gene_table, aes(x=Var1, y=Freq, fill = group))+
  geom_bar(stat="identity") +
  labs(x="Bolting regulator", y="Number of variants") +
  scale_fill_manual(values=c("orange", "red", "darkgrey")) + theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
snp_table <- as.data.frame(sort(table(sigres$SNP), decreasing = T)) %>%
  rowwise() %>%
  mutate(FDR = as.factor(ifelse(Var1 %in% maybelinked$SNP, maybelinked$FDR[which(maybelinked$SNP == Var1)], "other")))
p.snp_table <- ggplot(snp_table, aes(x=Var1, y=Freq))+
  geom_bar(stat="identity") +
  labs(x="Variant", y="Number of target genes") +
  scale_fill_manual(values=c("green", "darkgrey", 'purple',"blue", "red")) + theme_minimal() +
  theme(axis.text.x = element_blank())


ggsave(p.gene_table, file ="results/plots/gene_table.pdf", width= 140, height=90, units = "mm", dpi=300)
ggsave(p.snp_table, file ="results/plots/snp_table.pdf", width=190, height=90, units = "mm", dpi=300)
# pdf(file = "results/plots/expression_heatmap.pdf", width = 190 * 0.0393701, height = 110 * 0.0393701)
# print(expression_heatmap)
# dev.off()

# rm(list=ls())