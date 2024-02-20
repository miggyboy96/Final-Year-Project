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
  lm_model <- lm(data_long$expression ~ data_long$genotype)
  ggplot(data_long, aes(genotype, expression)) +
    geom_jitter(col="darkorange", alpha = 0.8, position=position_jitter(width=0.1)) +
    #geom_boxplot(outlier.size=0, alpha=0.6, fill="steelblue") +
    facet_grid(~comparison) +
    geom_smooth(method = 'lm', linewidth=1,col="darkred", aes(group=1), se=FALSE) +
    stat_cor(aes(group = 1,label =  paste(..p.label.., ..rr.label.., sep = "~~~~")))
}

calclog2FC <- function(snp, gene){
  # Extract expression values for the gene
  gene_expr <- expr[expr$geneid == gene, ]
  # Extract genotype values for the SNP
  variant_genotype <- gt[gt$snpid == snp, ]
  # Ensure matching samples between expression and genotype data
  common_samples <- intersect(names(gene_expr), names(variant_genotype))
  gene_expr <- gene_expr[common_samples]
  variant_genotype <- variant_genotype[common_samples]
   # Group expression values by genotype (0 vs 1 or 2)
  expr_genotype_before <- gene_expr[variant_genotype == min(variant_genotype)]
  expr_genotype_after <- gene_expr[variant_genotype == max(variant_genotype)] # either 1 or 2
  # Calculate average expression for each genotype group

  avg_expr_genotype_before <- mean(expr_genotype_before, na.rm = TRUE)
  avg_expr_genotype_after <- mean(expr_genotype_after, na.rm = TRUE)
  # Calculate log2 fold change
  log2fc <- log2(avg_expr_genotype_after / avg_expr_genotype_before)
  return(log2fc)
}

        
