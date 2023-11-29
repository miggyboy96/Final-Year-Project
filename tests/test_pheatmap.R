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
