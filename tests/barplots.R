##
#
#

## Packages
library(tidyr)
library(ggplot2)
library(ggpubr)

## Datasets
load(file = 'data/processed_data.Rda')
source(file =  'code/functions.R')

## Transposed versions of datasets
gt_trans <- data.frame(t(gt[,-1]))
colnames(gt_trans) <- t(gt[,1])
gt_trans <- tibble::rownames_to_column(gt_trans, 'sample')

expr_trans <- data.frame(t(expr[,-1]))
colnames(expr_trans) <- t(expr[,1])
expr_trans <- tibble::rownames_to_column(expr_trans, 'sample')

geno_long <- tidyr::gather(data = gt_trans, key = 'snp', value = 'genotype', -sample)
expr_long <- tidyr::gather(data = expr_trans, key = 'gene', value = 'expression', -sample)



snps <- names(which(snp_clusters==4))
genes <- names(which(gene_clusters==4))
snp.vs.gene(c(snps[rep(1:4,each=2)],snps[1]),c(rep(genes[c(2,4)],4),genes[1]))
snp.vs.gene(snps[rep(1:4,each=4)],genes[rep(1:4, 4)])
snp.vs.gene('snp1728', 'AT1G7880')