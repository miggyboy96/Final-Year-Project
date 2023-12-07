##
#
#

## Packages
library(tidyr)
library(ggplot2)

## Datasets
load(file = 'data/processed_data.Rda')

## Transposed versions of datasets
gt_trans <- data.frame(t(gt[,-1]))
colnames(gt_trans) <- t(gt[,1])
gt_trans <- tibble::rownames_to_column(gt_trans, 'sample')

expr_trans <- data.frame(t(expr[,-1]))
colnames(expr_trans) <- t(expr[,1])
expr_trans <- tibble::rownames_to_column(expr_trans, 'sample')

geno_long <- tidyr::gather(data = gt_trans, key = 'snp', value = 'genotype', -sample)
expr_long <- tidyr::gather(data = expr_trans, key = 'gene', value = 'expression', -sample)
data_long <- cbind(geno_long, expr_long["expression"])


snps <- gt$snpid[9:18]
genes <- expr$geneid[9:18]
