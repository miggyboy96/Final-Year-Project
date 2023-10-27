# Cleaning of data for eQTL

rm(list=ls())
wd <- "//userfs/ma2016/w2k/YEAR 3/Project/Data and Code"
setwd(wd)

# Online sources ####
# eQTL Tutorial https://rpubs.com/MajstorMaestro/349118
# Supplementary Materials https://www.biorxiv.org/content/10.1101/2023.09.11.557157v1.supplementary-material
# RNA-seq Data https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA1014093&o=acc_s%3Aa

# Packages ####

# Importing raw data and initial modifications ####
#Supplementary Materials were downloaded from the bioRxiv website into a subfolder of the working directory.

library(readr)

s11 <- read_tsv("Supplementary Materials/media-11.txt", skip = 8, show_col_types = F)
s12 <- readxl::read_excel("Supplementary Materials/media-12.xlsx", sheet = "TPM filtered", .name_repair = "minimal")

# Initial modifications #
gt <- as.data.frame(s11[-c(1,3:9)])
expr <- s12
colnames(expr)[1] <- "GENE"



# Genotype conversion to 0,1,2 ####

for (i in 1:nrow(gt)){
  for (j in 2:ncol(gt)){
    if (gt[i,j] == "0|0"){
      gt[i,j] = 0
    }
    else if (gt[i,j] == "1|1"){
      gt[i,j] = 2
    }
    else{
      gt[i,j] = 1
    }
  }
}

gt <- as.data.frame(sapply(gt, as.numeric))
gt <- gt[order(gt$"POS"),]



# Remove variants that are almost heterozygous across all samples. ####

index = c()
for (i in 1:nrow(gt)){
  if (length(grep(1, gt[i,-1]))>63){
    index <- append(index,i)}
}
gt <- gt[-index,]



# MAF threshold ####

gt.freq = rowMeans(gt[-1])/2 #allele (variant) freq.
threshold = 0.05 #the common threshold
gt <- gt[which(gt.freq > threshold),]



# SNP to gene dataframe ####

# Importing
tair10_gff <- read.table(url("https://arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"), sep="")
colnames(tair10_gff) <- c("Chr", "2", "region", "pos_start", "pos_end", "6","7","8","geneid")
tair10_gff <- subset(tair10_gff, region == "gene")

varinfo <- data.frame(pos=numeric(),geneid=character())
startrow = 1
for (snp in 1:length(gt[,1])){
  for (pos in length(tair10_gff[,1]):1){
    if (gt[snp,1] > tair10_gff[pos,"pos_start"]){
      newrow <- data.frame(pos=gt[snp,1], geneid=substring(tair10_gff[pos,"geneid"],4,12))
      varinfo <- rbind(varinfo,newrow)
      break
    }
  }
}
snpid = c()
for(variant in 1:length(varinfo[,1])){
  snpid <- c(snpid, paste("snp",variant))
}
 
varinfo <- cbind(snp,varinfo)     
gt[,1] <- snpid
colnames(gt)[1] <- "snpid"
colnames(expr)[1] <- "geneid"

# Saving variables ####

write.table(varinfo, file = "Data/snp_to_gene.txt", sep = "\t", row.names = F)
write.table(gt, file="Data/snp.txt", sep = "\t", row.names = F)
write.table(expr, file = "Data/expr.txt", sep = "\t", row.names = F)

save(list = c("gt", "expr", "varinfo", "wd"), file = "Data/eQTL_prep.rda")


