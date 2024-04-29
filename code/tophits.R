## tophits.R
## Miguel Alburo
## 28/04/2024
## Characterisation of top bolting variants

## Packages and data
library(tidyverse)
library(rvest)
load(file = "results/output_me.Rda")
load(file = "data/processed_data.Rda")

tophits <- distinct(sigres, SNP, .keep_all = T) %>%
  select(-pi) %>%
  rename(target = gene) %>%
  arrange(FDR) %>%
  mutate(gene = snp_to_gene[SNP]) %>%
  relocate(gene, .before = target)

### ... After obtaining gene descriptions and saving html result "Gene Description Search and Download" via TAIR (arabidopsis.org)

html.file <- "results/tophits_descriptions.html"
html.file %>%
  read_html %>%
  html_table() %>%
  .[[2]] -> tophits.descriptions

tophits.descriptions <- tophits.descriptions %>%
  rename(gene = X1, description = X3, alias = X5) %>%
  select(gene, description, alias) %>%
  slice(-1)

tophits <- inner_join(tophits,tophits.descriptions, by = "gene")
tophits <- tophits %>%
  relocate(SNP, gene, alias, description, target, FDR, log2FC) %>%
  arrange(FDR)

maybelinked <- tophits[duplicated(tophits$FDR) | duplicated(tophits$FDR, fromLast = TRUE),]

write.csv(tophits, file = "results/tophits.csv")