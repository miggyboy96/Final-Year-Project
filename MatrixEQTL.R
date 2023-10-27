
rm(list=ls())
wd <- "//userfs/ma2016/w2k/YEAR 3/Project/Data and Code"
setwd(wd)

load("Data/eQTL_prep.rda")

install.packages('IRkernel', repo='https://repo.miserver.it.umich.edu/cran/')
IRkernel::installspec()

#http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/
install.packages("MatrixEQTL")
library("MatrixEQTL")

# ####

covariates_file_name = character()
errorCovariance = numeric()

# ####

snps = SlicedData$new()
snps$fileDelimiter = "\t"
snps$fileOmitCharacters = "NA"
snps$fileSkipRows = 1
snps$fileSkipColumns = 1
snps$fileSliceSize = 2000
snps$LoadFile("Data/snp.txt")

expr = SlicedData$new()
expr$fileDelimiter = "\t"
expr$fileOmitCharacters = "NA"
expr$fileSkipRows = 1
expr$fileSkipColumns = 1
expr$fileSliceSize = 2000
expr$LoadFile("Data/expr.txt")

me = Matrix_eQTL_engine(
  snps = snps,
  gene = expr,
  output_file_name = "Data/MatrixEQTL_results.txt",
  pvOutputThreshold = 1e-2,
  useModel = modelLINEAR,
  verbose = F,
  errorCovariance = numeric()
)
