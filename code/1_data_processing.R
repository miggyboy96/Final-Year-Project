## 1_data_processing.R
# Miguel Alburo
# 25/11/2023
# Processing raw data from the supplementary materials into forms suitable for analyses.

## Load packages
library(tidyverse)
library(readxl)
source(file = 'code/functions.R')

### Importing raw data
## Expression and variant data. media-11.txt and media-12.txt downloaded from: https://www.biorxiv.org/content/10.1101/2023.09.11.557157v1.supplementary-material
S11 <- read_tsv(file = "data/raw/media-11.txt", skip = 8, show_col_types = F)
S12 <- readxl::read_excel(path = "data/raw/media-12.xlsx", sheet = "TPM filtered", .name_repair = "minimal")
GRN <- readxl::read_excel(path = 'data/raw/media-10.xlsx')

## Arabidopsis gene assembly data
tair10_gff <- read.table(url(description = "https://arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"), sep = "")

## Regulatory Genes
regulatorygenes <- unique(GRN$`regulatory gene`)

## Gene expression data
expr <- S12
colnames(expr)[1] <- "geneid"
geneids <- expr$geneid
sampleids <- colnames(expr)[-1]
write.table(expr,file = "data/expr.txt", sep = "\t", row.names = F)

## Arabidopsis genome annotation
geneloc <- tair10_gff %>% # Select and rename columns
  filter(V3 == "gene") %>%
  select(Chr = V1, left = V4, right = V5, geneid = V9)
geneloc$geneid <- sub(pattern = ".*Name=", replacement = "", geneloc$geneid) # Geneid substring from full string
write.table(geneloc,file = "data/geneloc.txt", sep = "\t", row.names = F)

# Manipulate chromosome data to use consistent naming
conversion_table <- data.frame(
  OriginalValue = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"),
  chr = c("1", "2", "3", "4", "5", "Pt", "Mt")
)
geneloc <- geneloc %>%
  left_join(conversion_table, by = c("Chr" = "OriginalValue"))
geneloc <- geneloc[c("geneid", "chr", "left", "right")]
geneloc$geneid <- as.character(geneloc$geneid)

## Variant genotype data
snpids <- paste0("snp_", seq_len(nrow(S11)))
snpsloc <- data.frame(
  snpid = snpids,
  chr = S11$`#CHROM`,
  pos = S11$POS
)
write.table(snpsloc,file = "data/snpsloc.txt", sep = "\t", row.names = F)
gt <- S11 %>%
  select(-one_of(c("#CHROM","FILTER", "INFO", "FORMAT" ,"POS", "QUAL", "ID", "REF", "ALT"))) %>%
  mutate(snpid = snpids, .before = 1)

# Conversion to 0,1,2 form
gt <- as.data.frame(gt)
for (i in seq_len(nrow(gt))){
  for (j in 2:ncol(gt)){
    if (gt[i,j] == "0|0"){
      gt[i,j] <- 0
    }
    else if (gt[i,j] == "1|1"){
      gt[i,j] <- 2
    }
    else{
      gt[i,j] <- 1
    }
  }
}
gt[,-1] <- as.data.frame(sapply(gt[,-1], as.numeric))
write.table(gt,file = "data/gt.txt", sep = "\t", row.names = F)

# Calculating minor allele frequency
maf <- sapply(snpids, function(snp){
  calculateMAF(gt[gt$snpid == snp, -1])
})
#min(maf) # = 0.0538

# Filtering SNPs based on occurances of 0,1,2
#threshold <- length(sampleids) - 2  # if heterozygous in all but 2 samples
#gt <- gt[-which(rowSums(gt[, -1] == 1) > threshold),]

## Creating a conversion table between snp and geneid
snp_to_gene <- mapply(pos.to.gene, snpsloc$pos, snpsloc$chr)
names(snp_to_gene) <- snpsloc$snpid

# Print results
print(paste(length(snpids),"variants/SNPs")) # no. variants
print(paste(length(sampleids),"samples")) # no. samples
print(paste(length(geneids),"genes")) # no. genes
print(paste(length(regulatorygenes),"bolting-associated regulatory genes")) # no. regulatory genes

## Saving and writing data

variables <- unique(c('gt', 'expr', 'regulatorygenes', 'snpids', 'geneids', 'sampleids', 'geneloc', 'snpsloc','snp_to_gene'))
save(list=variables, file = "data/processed_data.Rda")
rm(list=ls())