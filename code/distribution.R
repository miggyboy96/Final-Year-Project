## distribution.R
## Miguel Alburo
## 05/03/2024
##

# Packages and data
library(tidyverse)
library(rvest)
load(file = "results/output_me.Rda")
load(file = "data/processed_data.Rda")
source(file = "code/1_functions.R")


# Numeric vector of the end position of the previous chromosome. Used to create a single position metric.
boundary <- c( # The
  0,
  30425192,
  19696821 + 30425192,
  23459804 + 19696821 + 30425192,
  18584524 + 23459804 + 19696821 + 30425192,
  119136982
)
# chr_start_position <- function(x) {
#   if (all(is.numeric(x))) {
#     position_conversion <- data.frame(
#       chr = 1:5,
#       start = c(0,
#                 30425192,
#                 19696821 + 30425192,
#                 23459804 + 19696821 + 30425192,
#                 18584524 + 23459804 + 19696821 + 30425192)
#     )
#     # Use sapply to handle vector inputs
#     starts <- sapply(x, function(val) {
#       start_val <- position_conversion$start[position_conversion$chr == val]
#       if (length(start_val) == 0) return(NA)  # Return NA if no match is found
#       return(start_val)
#     })
#     return(starts)
#   } else {
#     stop("Input must be all integers.")
#   }
# }

# Converts chr,pos to single pos metric
snp_positions <- snpsloc %>%
  mutate(chr = as.integer(chr)) %>%
  filter(!is.na(chr)) %>%
  mutate(SNP_pos = pos + boundary[chr]) %>%
  select(snpid, SNP_pos) %>%
  rename(SNP = snpid)
gene_positions <- geneloc %>%
  mutate(chr = as.integer(chr)) %>%
  filter(!is.na(chr)) %>%
  mutate(gene_pos = (left + right )/ 2 + boundary[chr]) %>% ## midpoint of gene
  select(geneid, gene_pos) %>%
  rename(gene = geneid)

# long data of eQTL marker and target positions
locations_long <- sigres %>%
  mutate(log_FDR = -log(FDR,10)) %>%
  select(SNP, gene, log_FDR) %>%
  inner_join(snp_positions, by = "SNP") %>%
  inner_join(gene_positions, by = "gene") %>%
  select(SNP, gene, SNP_pos, gene_pos, log_FDR) %>%
  arrange(log_FDR)

# Cis trans plot prep
scale <- 1e-7 # Scaling position values
model <- lm(gene_pos ~ SNP_pos, data = locations_long)
r_squared <-  format(summary(model)$r.squared, scientific = T) # Extracting R-squared of cistrans plot

grid <- data.frame(
  xmin = rep(boundary[1:5], each = 5) * scale,
  xmax = rep(boundary[2:6], each = 5) * scale,
  ymin = rep(boundary[1:5], 5) * scale,
  ymax = rep(boundary[2:6], 5) * scale,
  fill = c(rep(c("blue","white"),12),"blue")
)



# Cis-trans plot
p.locations <- ggplot() +  # Map x and y to the dataframe columns
  geom_point(data = locations_long, aes(x = SNP_pos * scale, y = gene_pos * scale, colour = log_FDR), size = 0.8) + # Add points
  geom_rect(data = grid, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill)) +
  scale_fill_manual(values = alpha(c("black", "white"), 0.15), guide = "none") +
  geom_abline(slope = 1, intercept = 0, color = "blue", linetype = "dashed") +
  coord_cartesian(xlim = c(0, boundary[6] * scale), ylim = c(0, boundary[6] * scale), expand = FALSE)  +
  scale_x_continuous(breaks = boundary[-length(boundary)] + diff(boundary) / 2, # Midpoints of each chromosome
                     labels = 1:5) +
  scale_y_continuous(breaks = boundary[-length(boundary)] + diff(boundary) / 2, # Midpoints of each chromosome
                     labels = 1:5) +
  scale_color_gradient(low = "violet", high = "black", name = "log(FDR)") +
  annotate("text", x = Inf, y = Inf, label = paste0("R-squared = ", r_squared), hjust = 1.4, vjust = 3, size = 4, color = "red") +
  theme_bw() +
  labs(x = "SNP location", y = "Target regulator location")


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
tallydata$'Genome (Normalised)' <- table(loci2region(positionSNP(snpids))) * (mean(sapply(tallydata,mean)) / mean(table(loci2region(positionSNP(snpids))))) * 1.1

# Position of regulatory genes as 1MB regions
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

regionsdata <- read.table(file="data/genomic_regions.txt")

for (gene in seq_along(tallydata)){
  regions <- as.integer(names(tallydata[[gene]]))
  counts <- as.integer(tallydata[[gene]])
  tally_matrix[gene,regions] <- counts
}
tally_matrix <- as.data.frame(tally_matrix) %>%
  mutate(gene = rownames(.))

tally_long <- tally_matrix %>%
  pivot_longer(!gene, names_to = "region", values_to = "count") %>%
  mutate(region = as.integer(region)) %>%
  mutate(count = as.integer(count)) %>%
  rowwise() %>%
  mutate(chromosome = regionsdata$Chromosome[which(regionsdata$MbRegionLabel==region)])

# Add the chromosome column
tally_long$region <- as.integer(tally_long$region)
tally_long <- tally_long %>%
  rowwise() %>%
  mutate(chromosome = getChromosome(as.integer(region)))
tally_long$count <- as.integer(tally_long$count)

# Plot heatmap
p.eqtl_density <- ggplot(tally_long, aes(x = region, y = gene, fill = count)) +
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
    axis.text.x = element_text(angle = 90, hjust = 1),
  ) +
  coord_cartesian(xlim=c(0.5,121.5), expand = FALSE) +
  geom_vline(xintercept = c(31, 51, 75, 94)+0.5, linetype = "dashed", color = "black") +
  scale_x_continuous(breaks = c(16, 42, 64, 85, 108), # Midpoints of each region
                     labels = 1:5) +
  theme_minimal() +
  labs(x = "Chromosome", y = "Regulatory Gene", fill = "no. eQTL")  # Adjust axis labels as needed

html.file <- "results/at2g23320_araqtl.html" %>%
  read_html %>%
  html_table() %>%
  .[[1]] -> x
maybecorrelated <- as.vector(as.data.frame(sort(table(sigres$gene), decreasing = T))[1:10,1])
araqtl.res <- x[,-c(1,8)] %>%
  rename_at(1,~"tairid") %>%
    filter(tairid %in% maybecorrelated)



print(p.locations)
print(p.eqtl_density)
ggsave(p.locations, file = "results/plots/cistrans.pdf", width=120, height=90, units = "mm", dpi=300)
ggsave(p.eqtl_density, file = "results/plots/eqtl_hotspot.pdf", width=170, height=100, units = "mm", dpi=200)
save(list = c('tophits', 'locations_long', 'araqtl.res'), file = 'results/distribution.Rda')