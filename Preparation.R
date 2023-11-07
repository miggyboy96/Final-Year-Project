#### Preparing data for eQTL ####
# Clear the environment and set the working directory
rm(list = ls())
wd <- "//userfs/ma2016/w2k/YEAR 3/Project/Analysis"
setwd(wd)

# Load necessary libraries
library(readr)  # For reading data
library(dplyr)  # For data manipulation

#### Importing #####
#  Expression and variant data. media-11.txt and media-12.txt downloaded from: https://www.biorxiv.org/content/10.1101/2023.09.11.557157v1.supplementary-material
S11 <- read_tsv("Supplementary Materials/media-11.txt", skip = 8, show_col_types = F)
S12 <- readxl::read_excel("Supplementary Materials/media-12.xlsx", sheet = "TPM filtered", .name_repair = "minimal")

# Arabidopsis gene assembly data
tair10_gff <- read.table(url("https://arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"), sep = "")

#### Data manipulation ####
# Gene expression data
expr <- S12
colnames(expr)[1] <- "geneid"
geneid <- unlist(expr[,1],use.names = F)
sampleid <- colnames(expr)[-1]

# Arabidopsis gene assembly
geneloc <- tair10_gff %>%
  filter(V3 == "gene") %>%
  select(Chr = V1, left = V4, right = V5, geneid = V9)
geneloc$geneid <- sub(".*Name=", "", geneloc$geneid)

## Manipulate chromosome data to use consistent naming
conversion_table <- data.frame(
  OriginalValue = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "ChrC", "ChrM"),
  chr = c("1", "2", "3", "4", "5", "Pt", "Mt")
)
geneloc <- geneloc %>%
  left_join(conversion_table, by = c("Chr" = "OriginalValue"))

geneloc <- geneloc[c("geneid", "chr", "left", "right")]
geneloc$geneid <- as.character(geneloc$geneid)


# Variant data
snpid <- paste0("snp_", 1:nrow(S11))
snpsloc <- data.frame(
  snpid = snpid,
  chr = S11$`#CHROM`,
  pos = S11$POS
)
snp <- S11 %>% select(-one_of(c("#CHROM","FILTER", "INFO", "FORMAT" ,"POS", "QUAL", "ID", "REF", "ALT")))
snp <- cbind(snpid, snp)         

## Converting genotype data
for (i in 1:nrow(snp)){
  for (j in 2:ncol(snp)){
    if (snp[i,j] == "0|0"){
      snp[i,j] = 0
    }
    else if (snp[i,j] == "1|1"){
      snp[i,j] = 2
    }
    else{
      snp[i,j] = 1
    }
  }
}
snp[,-1] <- as.data.frame(sapply(snp[,-1], as.numeric))

## Removing snps heterozygous across nearly all samples
threshold = length(sampleid) - 2  # if heterozygous in all but 2 samples
snp <- snp[-which(rowSums(snp[, -1] == 1) > threshold),]

#### Saving and writing data ####
variables <- c('snp', 'expr', 'snpid', 'geneid', 'sampleid', 'geneloc', 'snpsloc')
save(list=variables, file = "Data/cleaned_data.rda")
for (x in variables){
  write.table(get(x),file = paste0("Data/",x,".txt"), sep = "\t", row.names = F)
}
