## conversion function to obtain geneid from genomic position
pos.to.gene <- function(pos, chr){
  gene <- NULL
  for (gene in 1:nrow(geneloc)){
    if (pos <= geneloc$right[gene] & chr == geneloc$chr[gene]){
      return(geneloc$geneid[gene])
      break
    }
  }
}
