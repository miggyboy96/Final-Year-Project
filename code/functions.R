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

snp.vs.gene <- function(snps,genes){
  data_long <- data.frame(sample = character(), snp = character(), genotype = numeric(), gene = character(), expression = numeric())
  for (index in seq(length(snps))){
    next_rows <- merge(subset(geno_long, snp == snps[[index]]), subset(expr_long, gene == genes[[index]]), by = 'sample')
    data_long <- rbind(data_long,next_rows)
  }
  data_long$comparison <- paste(data_long$snp,"vs.",data_long$gene)
  data_long$genotype <- as.factor(data_long$genotype)
  ggplot(data_long, aes(genotype, expression)) +
    geom_jitter(col="darkorange", position=position_jitter(width=0.25)) +
    geom_boxplot(outlier.size=0, alpha=0.6, fill="steelblue") +
    geom_smooth(method = 'lm', linewidth=1,col="darkred", aes(group=1), se=FALSE) +
    facet_wrap(~comparison) 
}
geom_sm



        
