## 2_matrixeqtl.R
## Miguel Alburo
## 25/11/2023
## Performing eQTL analysis using MatrixEQTL 

## Packages and data
library(tidyverse)
library("MatrixEQTL") # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
library(EnhancedVolcano)
library(ggsci)
source(file = "code/functions.R")
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

# Filtering based on regulatory genes
res <- output_eqtl %>%
  filter(gene %in% regulatorygenes) %>%
  select(SNP, gene, FDR) %>%  # Adjust these column names as per your data frame
  rowwise() %>%
  mutate(log2FC = calclog2FC(SNP,gene)) %>%
  mutate(pi = -1 * log10(FDR) * calclog2FC(SNP,gene)) %>%
  ungroup() %>%
  arrange(desc(abs(pi)))

## Filtering based on FDR and effect size
sigres <- res %>% filter(FDR<1e-5, (log2FC > 1 | log2FC < -1))
sigvariants <- unique(sigres$SNP)
sig_genes <- unique(sigres$gene)

# X-label
xlabel <- res$SNP
xlabel[duplicated(xlabel)] <- seq(1,sum(duplicated(xlabel)))

## Colours
#genetally <- count(sigres,gene,sort = T)
colour_key <- setNames(pal_ucscgb()(length(sig_genes)), sig_genes)
genecolours <- (colour_key[res$gene])
genecolours[is.na(genecolours)]  <- "gray30"
names(genecolours)[is.na(genecolours)] <- "NS"

# genecolours[which(! (res$gene %in% sig_genes))] <- "black"
# names(genecolours)[which(! (res$gene %in% sig_genes))] <- "p-value and log2FC"
# genecolours[which(res$FDR<1e-5 & (res$log2FC < 1 & res$log2FC > -1))] <- "gray30" #"royalblue"
# names(genecolours)[which(res$FDR<1e-5 & (res$log2FC < 1 & res$log2FC > -1))] <- "p-value"
# genecolours[which(res$FDR>1e-5 & (res$log2FC > 1 | res$log2FC < -1))] <- "gray30" # "forestgreen"
# names(genecolours)[which(res$FDR>1e-5 & (res$log2FC < 1 | res$log2FC > -1))] <- "log2FC"
# genecolours[which(res$FDR>1e-5 & (res$log2FC < 1 & res$log2FC > -1))] <- "gray30"
# names(genecolours)[which(res$FDR>1e-5 & (res$log2FC < 1 & res$log2FC > -1))] <- "NS"

# Draw Volcano plot
EnhancedVolcano(
  res,
  lab = xlabel,
  selectLab = unique(res$SNP[res$FDR < 1e-5]),
  x = 'log2FC',
  y = 'FDR',
  title = 'MatrixEQTL results',
  colCustom = genecolours,
  pCutoff = 1e-5,
  FCcutoff = 1,
  pointSize = 3.5,
  labSize = 6.0
)

# Print results
print(paste(sum(output_eqtl$FDR < 1e-5),"general eQTLs with FDR adjusted q<1e-5"))
print(paste(length(unique(output_eqtl$SNP[output_eqtl$FDR< 1e-5])), "unique variants are linked to"))
print(paste(length(unique(output_eqtl$gene[output_eqtl$FDR< 1e-5])), "unique regulatory genes affected"))

print(paste(nrow(res[res$FDR< 1e-5,]),"bolting-regulator linked significant eQTLs with q<1e-5 and log2FC >1/<-1"))
print(paste(length(sigvariants), "unique subset of significant variants are linked to"))
print(paste(length(sig_genes), "unique affected bolting-related regulatory genes"))

## Saving variables to .Rda
save(list = 'output_eqtl', 'res', "sigres", "sig_genes", "sigvariants", file = "results/output_me.Rda")
rm(list=ls())