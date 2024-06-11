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
  combined <- as.data.frame(t(full_join(genotype,expression, by = colnames(genotype))))
  colnames(combined) <- c("genotype", "expression")
  combined <- combined %>%
    arrange(genotype) %>%
    mutate(genotype = as.factor(as.integer(genotype)))
  return(as.data.frame(combined))
}

test <- genotype.expression("snp_42", "AT4G32800")

genotype.expression.plot <- function(snp, gene, bin_width){
  df <- genotype.expression(snp,gene)
  p <- ggplot(df, aes(x = expression, fill = genotype)) +
    geom_density(alpha = 0.4) +
    theme_bw() +
    coord_cartesian(expand = FALSE)  +
    labs(x = paste(gene, "expression (TPM)"), y = "density", fill = paste(snp, "genotype"))

  # for (g in df$genotype){
  #   hist_data <- filter(df,genotype == g) %>% select(expression) %>% unlist(as.vector(.)) %>% hist(x = ., plot = F)
  #   p <- p + ggplot()
  # }
  p
}