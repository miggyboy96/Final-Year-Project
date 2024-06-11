## go_enrichment.R
## Miguel Alburo
## 04/05/2024
## Gene ontology enrichment analysis on the genes of variants

## Load tidyverse and data
library(tidyverse)
source(file = "code/1_functions.R")
load(file = 'results/output_me.Rda')
load(file = 'data/processed_data.Rda')

# Specific packages
library(gprofiler2) # ORA

genelist <- snp_to_gene[unique(res$SNP)]

gsea <- gost(query = genelist, organism = "athaliana", ordered_query = TRUE)
p.gsea <- gostplot(gsea, interactive = F)
gost
# GO table and comparison with Ethan's results
ethansGO <- readxl::read_excel(path = 'data/raw/media-2.xlsx')
gotable <- as.data.frame(gsea$result) %>%
 select(-parents)



# Enrichment Plot and Table
publish_gostplot(p.gsea, filename = "results/plots/gsea.pdf", width = 0.0393701 * 140, height = 0.0393701 * 90)
publish_gosttable(gsea, ggplot = T, filename = "results/plots/gotable.pdf")
write.csv(gotable, file = "results/gotable.csv")