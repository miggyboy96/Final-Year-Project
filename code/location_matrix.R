## 6_snploc_heatmap.R
## Miguel Alburo
## 05/03/2024
##

# Packages
library(tidyverse)


# Load Data
load(file = "results/output_me.Rda")
load(file = "data/processed_data.Rda")
#load(file = "results/volcano.Rda")
source(file = "code/functions.R")

# gene | Mb region | frequency

# 1. For a given gene, returns a vector of significant SNP
# 2. For a vector of SNPs, returns genomic location and chromosome
# 3. For a data frame of chromosomes and genomic location, returns a vector of associated genomic region
# 4. For a vector of genomic regions, tallys frequency of each genomic region

## Constructing adjacency matrix
genelist <- unique(sigres$gene)
snplist <- unique(sigres$SNP)

## Single position metric conversion
chr_start_position <- function(x) {
  if (all(is.integer(x))) {
    position_conversion <- data.frame(
      chr = 1:5,
      start = c(0,30425192, 19696821, 23459804, 18584524)
    )
    # Use sapply to handle vector inputs
    starts <- sapply(x, function(val) {
      start_val <- position_conversion$start[position_conversion$chr == val]
      if (length(start_val) == 0) return(NA)  # Return NA if no match is found
      return(start_val)
    })
    return(starts)
  } else {
    stop("Input must be all integers.")
  }
}
snp_positions <- snpsloc %>%
  mutate(chr = as.integer(chr)) %>%
  filter(!is.na(chr)) %>%
  rowwise() %>%
  mutate(SNP_pos = pos + chr_start_position(chr)) %>%
  select(snpid, SNP_pos) %>%
  rename(SNP = snpid)
gene_positions <- geneloc %>%
  mutate(chr = as.integer(chr)) %>%
  filter(!is.na(chr)) %>%
  rowwise() %>%
  mutate(gene_pos = (left + right / 2) + chr_start_position(chr)) %>%
  select(geneid, gene_pos) %>%
  rename(gene = geneid)


# Location long data
locations_long <- sigres %>%
  mutate(log_FDR = -log(FDR,10)) %>%
  select(SNP, gene, log_FDR) %>%
  inner_join(snp_positions, by = "SNP") %>%
  inner_join(gene_positions, by = "gene") %>%
  select(SNP, gene, SNP_pos, gene_pos, log_FDR)

locations_plot <- ggplot(locations_long, aes(x = SNP_pos, y = gene_pos, color = log_FDR)) +  # Map x and y to the dataframe columns
  geom_point() +# Add points
  scale_color_gradient(low = "white", high = "red") +
  theme_bw() +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +# Use a minimal theme
  labs(x = "SNP location", y = "Target regulator location")       # Label the axes

print(locations_plot)
