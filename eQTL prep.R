# PREPARING DATA FOR eQTL

# Clearing environment and setting wd ####

rm(list=ls())
wd <- "//userfs/ma2016/w2k/YEAR 3/Project/Analysis"
setwd(wd)

# Importing raw data ####

# media-11.txt and media-12.txt downloaded from:
#https://www.biorxiv.org/content/10.1101/2023.09.11.557157v1.supplementary-material

library(readr) #tibble formats

S11 <- read_tsv("Supplementary Materials/media-11.txt", skip = 8, show_col_types = F)
S12 <- readxl::read_excel("Supplementary Materials/media-12.xlsx", sheet = "TPM filtered", .name_repair = "minimal")
tair10_gff <- read.table(url("https://arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"), sep="")

# Initial modifications ####

# Gene expression data
expr <- S12
colnames(expr)[1] <- "geneid"
write.table(expr,file="Data/expr.txt",sep="\t")

# Arabidopsis gene assembly
tair10 <- subset(tair10_gff, V3 == "gene")
tair10 <- tair10[,-c(2,3,(6:8))]
colnames(tair10) <- c("Chr", "pos_start", "pos_end", "geneid")
tair10$geneid <- lapply(tair10$geneid, function (x) {sub(".*Name=","",x)})

# SNP information
snps <- data.frame(snpid=character(), Chr=character(), pos=numeric(), geneid=character())
for (snp in 1:nrow(S11)){
  nextpos = S11$POS[snp]
  nextsnp = paste("snp_",snp,sep="")
  nextChr = S11$`#CHROM`[snp]
  for (pos in nrow(tair10):1){
    if (nextpos > tair10$pos_start[pos]){
      nextgeneid = tair10$geneid[pos]
      break
    }
  }
  nextrow = data.frame(snpid=nextsnp,Chr=nextChr,pos=nextpos,geneid=nextgeneid)
  snps = rbind(snps,nextrow)
}

gtyp <- as.data.frame(s11[-c(1,3:9)]) #removing unneeded columns
