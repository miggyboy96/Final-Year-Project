## 2_matrixeqtl.R
## Miguel Alburo
## 25/11/2023
## Performing eQTL analysis using MatrixEQTL 

## Packages and data
library(tidyverse)
source(file = "code/1_functions.R")
library("MatrixEQTL") # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
library(EnhancedVolcano)
library(ggsci)
load(file = "data/processed_data.Rda")

## Preparing for analysis
# Linear model to use. Either modelLINEAR, modelANOVA or modelLINEAR_CROSS
useModel <- modelLINEAR
# Genotype file name
SNP_file_name <- "data/gt.txt"
snps_location_file_name <- "data/snpsloc.txt"
# Gene expression file name
expression_file_name <- "data/expr.txt"
gene_location_file_name <- "data/geneloc.txt"
# Covariates file name. Set to character() for no covariates
covariates_file_name <- character()
# Output file name
#output_file_name_cis <- "results/cis-output_me.txt"
output_file_name <- "results/output_me.txt"
# Only associations significant at this level will be saved
pvOutputThreshold_cis <- 2e-2
pvOutputThreshold <- 1e-2
# Error covariance matrix. Set to numeric() for identity.
errorCovariance <- numeric()
# Distance for local gene-SNP pairs
#cisDist <- 1e6 # Wihin 1 Megabase
## Load genotype data
snps <- SlicedData$new()
snps$fileDelimiter <- "\t" # the TAB character
snps$fileOmitCharacters <- "NA" # denote missing values
snps$fileSkipRows <- 1 # one row of column labels
snps$fileSkipColumns <- 1 # one column of row labels
snps$fileSliceSize <- 2000 # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)
## Load gene expression data
gene <- SlicedData$new()
gene$fileDelimiter <- "\t" # the TAB character
gene$fileOmitCharacters <- "NA" # denote missing values
gene$fileSkipRows <- 1 # one row of column labels
gene$fileSkipColumns <- 1 # one column of row labels
gene$fileSliceSize <- 2000 # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)
## Load covariates
cvrt <- SlicedData$new()

## Run analysis
## Import snpspos genepos
snpspos <- read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos <- read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
me <- Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  output_file_name = output_file_name,
  pvOutputThreshold = pvOutputThreshold,
  useModel = useModel,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  snpspos = snpspos,
  genepos = genepos,
  pvalue.hist = TRUE,
  min.pv.by.genesnp = FALSE,
  noFDRsaveMemory = FALSE
  #cisDist = cisDist,
  #output_file_name.cis = output_file_name_cis,
  #pvOutputThreshold.cis = pvOutputThreshold_cis,
)
# unlink(output_file_name_tra)
# unlink(output_file_name_cis)

## Reading Matrix_EQTL results
#output_eqtl_cis <- read.table(file = "results/cis-output_me.txt", header = T)
#output_eqtl_trans <- read.table(file = "results/trans-output_me.txt", header = T)
output_eqtl <- read.table(file = "results/output_me.txt", header = T) #rbind(output_eqtl_cis,output_eqtl_trans)

df <- output_eqtl %>%
  select(SNP, gene, FDR) %>%
  arrange(FDR) %>%
  filter(!duplicated(SNP)) %>%
  rowwise() %>%
  mutate(log2FC = calclog2FC(SNP,gene)) %>%
  ungroup()

p.volcano.all <- EnhancedVolcano(
  df,
  lab = df$SNP,
  selectLab = df$SNP[1:10],
  x = 'log2FC',
  xlim = c(-10,10),
  ylim = c(0,40),
  y = 'FDR',
  title = 'MatrixEQTL results',
  legendLabSize = 8,
  legendIconSize = 3.0,
  pCutoff = 1e-5,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0
)
print(p.volcano.all)
ggsave(p.volcano.all, file ="results/plots/volcanoplot_all.pdf", width=190, height=150, units = "mm", dpi=300)
# Print results
print(paste(sum(output_eqtl$FDR < 1e-5),"general eQTLs with FDR adjusted q<1e-5"))
print(paste(length(unique(output_eqtl$SNP[output_eqtl$FDR< 1e-5])), "unique variants are linked to"))
print(paste(length(unique(output_eqtl$gene[output_eqtl$FDR< 1e-5])), "unique genes affected"))


rm(list=ls())