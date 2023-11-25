## eqtl_analysis.R
## Miguel Alburo
## 25/11/2023
## Performing eQTL analysis using MatrixEQTL 

## Load packages
  library("MatrixEQTL") # http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/

## Import processed data
  load("data/processed_data.Rda")
  
## Preparing for analysis
  # Linear model to use. Either modelLINEAR, modelANOVA or modelLINEAR_CROSS
    useModel = modelLINEAR
  # Genotype file name
    SNP_file_name = "data/gt.txt";
    snps_location_file_name = "data/snpsloc.txt";
  # Gene expression file name
    expression_file_name = "data/expr.txt";
    gene_location_file_name = "data/geneloc.txt";
  # Covariates file name. Set to character() for no covariates
    covariates_file_name = character();
  # Output file name
    output_file_name_cis = "data/cis-output_me.txt";
    output_file_name_tra = "data/trans-output_me.txt";
  # Only associations significant at this level will be saved
    pvOutputThreshold_cis = 2e-2;
    pvOutputThreshold_tra = 1e-2;
  # Error covariance matrix. Set to numeric() for identity.
    errorCovariance = numeric();
  # Distance for local gene-SNP pairs
    cisDist = 1e6;
  ## Load genotype data
    snps = SlicedData$new();
    snps$fileDelimiter = "\t"; # the TAB character
    snps$fileOmitCharacters = "NA"; # denote missing values;
    snps$fileSkipRows = 1; # one row of column labels
    snps$fileSkipColumns = 1; # one column of row labels
    snps$fileSliceSize = 2000; # read file in slices of 2,000 rows
    snps$LoadFile(SNP_file_name);
  ## Load gene expression data
    gene = SlicedData$new();
    gene$fileDelimiter = "\t"; # the TAB character
    gene$fileOmitCharacters = "NA"; # denote missing values;
    gene$fileSkipRows = 1; # one row of column labels
    gene$fileSkipColumns = 1; # one column of row labels
    gene$fileSliceSize = 2000; # read file in slices of 2,000 rows
    gene$LoadFile(expression_file_name);
  ## Load covariates
    cvrt = SlicedData$new()
    
## Run analysis
  ## Import snpspos genepos
    snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE);
    genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE);
  me = Matrix_eQTL_main(
    snps = snps,
    gene = gene,
    cvrt = cvrt,
    output_file_name = output_file_name_tra,
    pvOutputThreshold = pvOutputThreshold_tra,
    useModel = useModel,
    errorCovariance = errorCovariance,
    verbose = TRUE,
    output_file_name.cis = output_file_name_cis,
    pvOutputThreshold.cis = pvOutputThreshold_cis,
    snpspos = snpspos,
    genepos = genepos,
    cisDist = cisDist,
    pvalue.hist = TRUE,
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE
  )
  # unlink(output_file_name_tra);
  # unlink(output_file_name_cis);


## Saving variables to .Rda ##
me.cis <- read.table("data/cis-output_me.txt", header=T)
me.trans <- read.table("data/trans-output_me.txt", header=T)
save(list=c('me.cis', 'me.trans'), file="data/output_me.Rda")
