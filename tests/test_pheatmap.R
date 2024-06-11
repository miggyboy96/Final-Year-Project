load('data/processed_data.Rda')
load('../results/output_me.Rda')



## Loading Ethan's GRN
  library(readxl)
  grn <- readxl::read_excel('data/raw/media-10.xlsx')
  
grn.genes <- unique(c(grn$`regulatory gene`,grn$`target gene`))
regulatoryGenes <- unique(grn$`regulatory gene`)

meLong <- rbind(me.cis,me.trans)

library(dplyr)

subsetmeLong <- meLong[which(meLong$gene %in% regulatoryGenes),]
connectivity <- as.data.frame(table(subsetmeLong$SNP))
subsub <- as.vector(connectivity[which(connectivity$Freq>38),]$Var1)
subsetmeLong <- subsetmeLong[which(subsetmeLong$SNP %in% subsub), ]

library(pheatmap)
library(RColorBrewer)
library(viridis)


adjacencyMatrix <- matrix(0,
                          nrow = length(unique(subsetmeLong$SNP)),
                          ncol = length(unique(subsetmeLong$gene)),
                          dimnames = list(unique(subsetmeLong$SNP), unique(subsetmeLong$gene))
)

for (edge in seq_len(nrow(subsetmeLong))){
  snp <- subsetmeLong$SNP[edge]
  gene <- subsetmeLong$gene[edge]
  logp <- -log(subsetmeLong$p.value[edge], base=10)
  adjacencyMatrix[snp,gene] <- logp
}

pheatmap(t(adjacencyMatrix), color = inferno(60))


library(igraph)
GRN <- readxl::read_excel(path = 'data/raw/media-10.xlsx')
p.grn <- graph_from_data_frame(GRN, directed = TRUE)
plot(p.grn,
     vertex.size = 20,               # Controls the size of the nodes
     vertex.label = NULL,  # Node labels
     vertex.label.cex = 0.8,         # Font size for node labels
     vertex.color = "skyblue",       # Node color
     edge.arrow.size = 0.5,          # Arrow size
     edge.curved = 0.1,              # Edge curvature
     main = "Gene Regulatory Network")  # Title of the plot

ethansGO <- readxl::read_excel(path = 'data/raw/media-2.xlsx') %>%
  select(source, term_name, term_id, term_size, p_values_up, p_values_down) %>%
  group_by(source)

bplist <- ethansGO