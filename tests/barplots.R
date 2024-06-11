##
#
#

## Packages
library(tidyverse)

## Datasets
load(file = 'data/processed_data.Rda')
source(file = 'code/1_functions.R')


genotype.expression <- function(snp, gene){
  genotype <- filter(gt, snpid == snp) %>%
    select(-snpid)
  expression <- filter(expr, geneid == gene) %>%
    select(-geneid)
  combined <- t(full_join(genotype,expression, by = colnames(genotype)))
  colnames(combined) <- c("genotype", "expression")
  return(combined)}
## Transposed versions of datasets
gt_trans <- data.frame(t(gt[,-1]))
colnames(gt_trans) <- t(gt[,1])
gt_trans <- tibble::rownames_to_column(gt_trans, 'sample')

expr_trans <- data.frame(t(expr[,-1]))
colnames(expr_trans) <- t(expr[,1])
expr_trans <- tibble::rownames_to_column(expr_trans, 'sample')

geno_long <- tidyr::gather(data = gt_trans, key = 'snp', value = 'genotype', -sample)
expr_long <- tidyr::gather(data = expr_trans, key = 'gene', value = 'expression', -sample)



snps <- c("snp_42", "snp_1841")
genes <- sig_genes[1:5]
snp.vs.gene(snps[rep(1:2,each=5)],genes[rep(1:5,2)])
snp.vs.gene(snps[rep(1:4,each=4)],genes[rep(1:4, 4)])
snp.vs.gene('snp1728', 'AT1G7880')

snp_42_log2FC<- unlist(lapply(sig_genes, function(x){calclog2FC("snp_42",x)}))
snp_1841_log2FC<- unlist(lapply(sig_genes, function(x){calclog2FC("snp_1841",x)}))
df <- data.frame(x = seq(length(sig_genes)), snp_42 = snp_42_log2FC, snp_1841 = snp_1841_log2FC)

ggplot(df, aes(x)) +                    # basic graphical object
  geom_line(aes(y=snp_42), colour="red") +  # first layer
  geom_line(aes(y=-snp_1841), colour="green")  # second layer
plot(1:20, snp_42_log2FC)
plot(1:20, -snp_1841_log2FC)
