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

tallydata <- sapply(genelist, function(x){
  relatedSNPs <- relatedSNP(x)
  positions <- positionSNP(relatedSNPs)
  SNPregions <- loci2region(positions)
  tallys <- table(SNPregions)
  return(tallys)
})

control_genome <- table(loci2region(positionSNP(snpids)))
control_genome <- control_genome * (mean(sapply(tallydata,mean)) / mean(control_genome)) * 1.1 # Rescale
tallydata$'Genome (Normalised)' <- control_genome

# Position of regulatory genes
genepos <- geneloc %>%
  select(-right) %>%
  filter(geneid %in% genelist) %>%
  rename(pos = left)
generegions <- genepos %>%
  mutate(region = loci2region(genepos)) %>%
  select(geneid, region) %>%
  rename(gene = geneid)

tally_matrix <- matrix(
  data = 0,
  dimnames = list(names(tallydata), regionsdata$MbRegionLabel),
  nrow = length(tallydata), ncol = 121
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
  geom_tile(data = generegions, color = "red", fill = NA) +
  scale_fill_gradient(low = "white", high = "blue", breaks = seq(10,0,-5), limits = c(0,10), guide = "colorbar") +  # Adjust colors as needed
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


### location matrix plot

# location matrix
location_matix <- output_eqtl %>%
  mutate(log_fdr = - log(FDR, 10)) %>%
  select(SNP, gene, log_fdr) %>%

# position of ALL genes
genepos <- geneloc %>%
  select(-right) %>%
  filter(geneid %in% genelist) %>%
  rename(pos = left)
generegions <- genepos %>%
  mutate(region = loci2region(genepos)) %>%
  select(geneid, region) %>%
  rename(gene = geneid)

