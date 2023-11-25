## data_processing.R
## Miguel Alburo
## 25/11/2023
## Processing raw data from the supplementary materials into forms suitable for analyses.

## Load packages
  library(readr)  # For reading data
  library(readxl)
  library(dplyr)  # For data manipulation

## Importing raw data
  ## Expression and variant data. media-11.txt and media-12.txt downloaded from: https://www.biorxiv.org/content/10.1101/2023.09.11.557157v1.supplementary-material
    S11 <- read_tsv("data/raw/media-11.txt", skip = 8, show_col_types = F)
    S12 <- readxl::read_excel("data/raw/media-12.xlsx", sheet = "TPM filtered", .name_repair = "minimal")
  ## Arabidopsis gene assembly data
    tair10_gff <- read.table(url("https://arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff"), sep = "")

## Gene expression data
  expr <- S12
  colnames(expr)[1] <- "geneid"
  geneid <- unlist(expr[,1],use.names = F)
  sampleid <- colnames(expr)[-1]

## Arabidopsis genome annotation
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

## Variant genotype data
  snpid <- paste0("snp_", 1:nrow(S11))
  snpsloc <- data.frame(
    snpid = snpid,
    chr = S11$`#CHROM`,
    pos = S11$POS
  )
  gt <- S11 %>% select(-one_of(c("#CHROM","FILTER", "INFO", "FORMAT" ,"POS", "QUAL", "ID", "REF", "ALT")))
  gt <- cbind(snpid, gt)         
  ## Conversion to 0,1,2 form
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
    gt[,-1] <- as.data.frame(sapply(gt[,-1], as.numeric))
  ## Removing snps heterozygous across nearly all samples
    threshold = length(sampleid) - 2  # if heterozygous in all but 2 samples
    gt <- gt[-which(rowSums(gt[, -1] == 1) > threshold),]

## Saving and writing data
  variables <- unique(c('variables','gt', 'expr', 'snpid', 'geneid', 'sampleid', 'geneloc', 'snpsloc'))
  save(list=variables, file = "data/processed_data.Rda")
  for (x in variables){
    write.table(get(x),file = paste0("data/",x,".txt"), sep = "\t", row.names = F)
  }