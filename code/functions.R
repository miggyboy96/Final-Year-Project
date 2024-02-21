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

myManhattan <- function(df, graph.title = "", highlight = NULL, highlight.col = "green",
                        col = c("lightblue", "navy"), even.facet = FALSE, chrom.lab = NULL,
                        suggestiveline = 1e-05, suggestivecolor = "blue",
                        genomewideline = 5e-08, genomewidecolor = "red",
                        font.size = 12, axis.size = 0.5, significance = NULL, report = FALSE,
                        inf.corr = 0.95, y.step = 2, point.size = 1, label.snps = NULL){
  myMin <- min(df$P[df$P != 0]) * inf.corr
  df$P[df$P == 0] <- myMin
  require(ggplot2)
  require(ggrepel)
  require(stats)
  y.title <- expression(-log[10](italic(p)))
  if (length(col) > length(unique(df$CHR))){
    chrom.col <- col[1:length(unique(df$CHR))]
  } else if (!(length(col) > length(unique(df$CHR)))){
    chrom.col <- rep(col, length(unique(df$CHR))/length(col))
    if (length(chrom.col) < length(unique(df$CHR))){
      dif <- length(unique(df$CHR)) - length(chrom.col)
      chrom.col <- c(chrom.col, col[1:dif])
    }
  }
  y.max <- floor(max(-log10(df$P))) + 1
  if (y.max %% 2 != 0){
    y.max <- y.max + 1
  }
  if (!is.null(chrom.lab)){
    if (length(unique(df$CHR)) != length(chrom.lab)){
      warning("Number of chrom.lab different of number of chromosomes in dataset, argument ignored.")
    } else {
      df$CHR <- factor(df$CHR, levels = unique(df$CHR), labels=chrom.lab)
    }
  }
  g <- ggplot(df) +
    geom_point(aes(BP, -log10(P), colour = as.factor(CHR)), size = point.size)
  if (!is.null(significance)){
    if (is.numeric(significance)){
      genomewideline <- significance
      suggestiveline <- genomewideline / 0.005
    } else if (significance == "Bonferroni"){
      BFlevel <- 0.05 / length(df$SNP)
      cat("Bonferroni correction significance level:", BFlevel, "\n")
      genomewideline <- BFlevel
      suggestiveline <- BFlevel / 0.005
    } else if (significance == "FDR"){
      df$fdr <- p.adjust(df$P, "fdr")
      genomewideline <- 0.05
      suggestiveline <- FALSE
      y.title <- expression(-log[10](italic(q)))
      g <- ggplot(df) +
        geom_point(aes(BP, -log10(fdr), colour = as.factor(CHR)), size = point.size)
      if (!is.null(highlight)) {
        if (is.numeric(highlight)){
          highlight <- as.character(df$SNP[df$P < highlight])
        }
        if (any(!(highlight %in% df$SNP))){
          warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
        } else {
          g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                              aes(BP, -log10(fdr), group=SNP, colour=SNP),
                              color = highlight.col, size = point.size)
          highlight <- NULL
          y.max <- floor(max(-log10(df$fdr))) + 1
          if (y.max %% 2 != 0){
            y.max <- y.max + 1
          }
        }
      }
    }
  }
  if (even.facet){
    g <- g + facet_grid(.~CHR, scale = "free_x", switch = "x")
  } else {
    g <- g + facet_grid(.~CHR, scale = "free_x", space = "free_x", switch = "x")
  }
  g <- g + scale_colour_manual(values = chrom.col) +
    scale_y_continuous(expand = c(0, 0), limit = c(0, y.max),
                       breaks = seq(from = 0, to = y.max, by = y.step)) +
    scale_x_continuous() +
    theme(strip.background = element_blank(), legend.position = "none",
          axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.spacing.x=unit(0.1, "lines"),
          axis.line.y = element_line(size = axis.size, color = "black"),
          axis.ticks.y = element_line(size = axis.size, color = "black"),
          axis.ticks.length = unit(axis.size * 10, "points"),
          plot.title = element_text(hjust = (0.5), size = font.size + 8),
          axis.title.y = element_text(size = font.size + 5),
          axis.title.x = element_text(size = font.size + 5),
          axis.text = element_text(size = font.size),
          strip.text.x = element_text(size = font.size))+
    labs(title = graph.title, x = "Chromosome", y = y.title)
  if (!is.null(highlight)) {
    if (is.numeric(highlight)){
      highlight <- as.character(df$SNP[df$P < highlight])
    }
    if (any(!(highlight %in% df$SNP))){
      warning("Cannot highlight SNPs not present in the dataset. Argument is ignored.")
    } else {
      g <- g + geom_point(data = df[which(df$SNP %in% highlight), ],
                          aes(BP, -log10(P), group=SNP, colour=SNP),
                          color = highlight.col, size = point.size)
    }
  }
  if (suggestiveline){
    g <- g + geom_hline(yintercept = -log10(suggestiveline), color = suggestivecolor)
  }
  if (genomewideline){
    g <- g + geom_hline(yintercept = -log10(genomewideline), color = genomewidecolor)
  }
  if (report){
    if (significance == "FDR"){
      rep <- df[df$fdr < 0.05, ]
    } else if (significance == "Bonferroni"){
      rep <- df[df$P < BFlevel, ]
    } else if (is.numeric(significance)){
      rep <- df[df$P < significance, ]
    } else {
      cat("using default significance level, 5e-8")
      rep <- df[df$P < 5e-8, ]
    }
    print(rep)
  }
  # Load the ggrepel package
  require(ggrepel)

  # Check if label.snps is not null and is a vector
  if (!is.null(label.snps) && is.vector(label.snps)) {
    # Filter the data frame to only include SNPs that need to be labelled
    label.data <- df[df$SNP %in% label.snps,]

    # Add labels to the plot using ggrepel
    g <- g + geom_text_repel(data = label.data,
                             aes(x = BP, y = -log10(P), label = SNP),
                              font.size = 5)
#, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50')
  }

  return(g)
}

        
