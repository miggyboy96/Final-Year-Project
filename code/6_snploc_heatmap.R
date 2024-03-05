## 6_snploc_heatmap.R
## Miguel Alburo
## 05/03/2024
##

# Packages
library(tidyverse)

# gene | Mb region | frequency

# 1. For a given gene, returns a vector of significant SNP
# 2. For a vector of SNPs, returns genomic location and chromosome
# 3. For a data frame of chromosomes and genomic location, returns a vector of associated genomic region
# 4. For a vector of genomic regions, tallys frequency of each genomic region

# 1.
relatedSNP <- function(genes) {
  sapply(genes, function(x){
    output_eqtl[output_eqtl$gene==x, 1]
  })
}

# 2.
positionSNP <- function(snps) {
  positions <- snpsloc[snpsloc$snpid %in% snps,]
}

# 3.
regionsdata <- read.table(file="data/genomic_regions.txt")
loci2region <- function(loci_df) {
  # Loop over each locus in loci_df and find the corresponding MbRegion from chromosome_regions
  region <- as.integer(sapply(1:nrow(loci_df), function(idx) {
    # Extract the current locus chromosome and position
    current_chr <- loci_df$chr[idx]
    current_pos <- loci_df$pos[idx]

    # Find the corresponding MbRegion
    mb_region <- regionsdata$MbRegionLabel[regionsdata$Chromosome == current_chr &
                                                  regionsdata$RegionStart <= current_pos &
                                                  (regionsdata$RegionStart + 1e6 - 1) >= current_pos]
    return(mb_region)
  }))
  return(region)
}

## Constructing adjacency matrix
genelist <- unique(sigres$gene)

tallydata <- sapply(genelist, function(x){
  relatedSNPs <- relatedSNP(x)
  positions <- positionSNP(relatedSNPs)
  SNPregions <- loci2region(positions)
  tallys <- table(SNPregions)
  return(tallys)
})

tally_matrix <- matrix(
  data = 0,
  dimnames = list(genelist, regionsdata$MbRegionLabel),
  nrow = length(genelist), ncol = 121
)

for (gene in seq_along(tallydata)){
  regions <- as.integer(names(tallydata[[gene]]))
  counts <- as.integer(tallydata[[gene]])
  tally_matrix[gene,regions] <- counts
}
tally_matrix <- as.data.frame(tally_matrix)
tally_matrix$gene <- rownames(tally_matrix)

tally_long <- tally_matrix %>%
  pivot_longer(!gene, names_to = "region", values_to = "count")

# Add the chromosome column
tally_long$region <- as.integer(tally_long$region)
tally_long <- tally_long %>%
  rowwise() %>%
  mutate(chromosome = getChromosome(as.integer(region)))
tally_long$count <- as.integer(tally_long$count)

# Plot heatmap

ggplot(tally_long, aes(x = region, y = gene, fill = count)) +
  geom_tile(color = "gray90", alpha = 1) +  # This creates the heatmap tiles
  scale_fill_gradient(low = "white", high = "blue") +  # Adjust colors as needed
  theme(
    #bold font for legend text
    legend.text=element_text(face="bold"),
    #set thickness of axis ticks
    axis.ticks=element_line(size=0.4),
    #remove plot background
    plot.background=element_blank(),
    #remove plot border
    panel.border=element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1)
  ) +
  xlim(0.5,121.5) +
  geom_vline(xintercept = c(31, 51, 75, 94)+0.5, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(16, 42, 64, 85, 108), # Midpoints of each region
                     labels = 1:5) +
  theme_minimal() +
  labs(x = "Chromosome", y = "Regulatory Gene", fill = "SNP frequency")  # Adjust axis labels as needed
adj_matrix_long <- tally_long %>%
  mutate(chr = cut(as.numeric(region), breaks = seq(0, 121, by = 10),
                     labels = paste(seq(1, 12, by = 1),
                                    seq(10, 120, by = 10), sep = "-"),
                     include.lowest = TRUE))