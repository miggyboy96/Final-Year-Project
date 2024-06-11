## 3_tophits.R
## Miguel Alburo
## 28/04/2024
## Characterisation of top bolting variants

## Packages and data
library(tidyverse)
library(EnhancedVolcano)
library(rvest)
library(ggsci)
library(gridExtra)
load(file = 'data/processed_data.Rda')
source(file = 'code/1_functions.R')

## load eQTL results
output_eqtl <- read.table(file = "results/output_me.txt", header = T)

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
genecolours[which(res$FDR>1e-5 | (res$log2FC < 1 & res$log2FC > -1))] <- "gray30"
names(genecolours)[which(res$FDR>1e-5 | (res$log2FC < 1 & res$log2FC > -1))] <- "NS"

# Volcano Plot
# Draw Volcano plot
p.volcano <- EnhancedVolcano(
  res,
  lab = xlabel,
  selectLab = unique(res$SNP[res$FDR < 1e-5]),
  x = 'log2FC',
  y = 'FDR',
  title = 'MatrixEQTL results',
  colCustom = genecolours,
  legendLabSize = 8,
  legendIconSize = 3.0,
  legendPosition = "none",
  pCutoff = 1e-5,
  FCcutoff = 1,
  pointSize = 2.0,
  labSize = 4.0
)
# Detailed significant results
tophits <- distinct(sigres, SNP, .keep_all = T) %>%
  select(-pi) %>%
  rename(target = gene) %>%
  arrange(FDR) %>%
  mutate(gene = snp_to_gene[SNP]) %>%
  relocate(gene, .before = target)

### ... After obtaining gene descriptions and saving html result "Gene Description Search and Download" via TAIR (arabidopsis.org)

html.file <- "results/tophits_descriptions.html"
html.file %>%
  read_html %>%
  html_table() %>%
  .[[2]] -> tophits.descriptions

tophits.descriptions <- tophits.descriptions %>%
  rename(gene = X1, description = X3, alias = X5) %>%
  select(gene, description, alias) %>%
  slice(-1)

tophits <- inner_join(tophits,tophits.descriptions, by = "gene")
tophits <- tophits %>%
  relocate(SNP, gene, alias, description, target, FDR, log2FC) %>%
  arrange(FDR)
write.csv(tophits, file = "results/tophits.csv")

df <- tophits %>%
  select(-description) %>%
  slice(1:20)
colours <- matrix("white", nrow(df), ncol(tophits[1:20,-4]))
fontcolours <- matrix("black", nrow(df), ncol(df))
fontfaces <- matrix("plain", nrow(df), ncol(df))
th <- gridExtra::ttheme_default(base_size = 10,
                                padding = grid::unit(c(4, 4), "mm"),
                                core=list(
                                  padding.h = grid::unit(c(15,15), "mm"),
                                  padding.v = grid::unit(c(15,15), "mm"),
                                  bg_params = list(fill = colours, col="black", lwd = 0.5),
                                  fg_params=list(hjust = 0, x = 0.02, col=fontcolours, fontface=fontfaces)),
                                colhead=list(bg_params = list(fill = "gray99", lwd = 0.5, col = "black"),
                                             fg_params=list(col="gray39", fontface="bold")),
                                rowhead=list(fg_params=list(col="black", fontface="bold")))
tb <- tableGrob(df, theme = th, rows = NULL)
h <- grid::unit.c(sum(tb$heights))
w <- grid::unit.c(sum(tb$widths))
tg <- gridExtra::arrangeGrob(tb, ncol = 1, widths = w, heights = h, bottom = grid::textGrob("tophits", x = 0.95, hjust = 1, gp = grid::gpar(fontsize=10, font=8, col = "cornflowerblue")))
width <- grid::convertWidth(sum(tg$widths), "in", TRUE) + 0.2
height <- grid::convertHeight(sum(tg$heights), "in", TRUE) + 0.2
p.tg <- ggplot2::ggplot() + ggplot2::annotation_custom(tg) + ggplot2::geom_blank() + ggplot2::theme_void()
ggplot2::ggsave(filename = "results/plots/tophits.pdf", plot = p.tg, height = height, width = width, limitsize = F)


# Linked variants
maybelinked <- sigres[duplicated(sigres$FDR) | duplicated(sigres$FDR, fromLast = TRUE),] %>%
  distinct(SNP, .keep_all = T)

# Saving variables
print(p.volcano)
ggsave(p.volcano, file ="results/plots/volcanoplot.pdf", width=190, height=150, units = "mm", dpi=300)
print(paste(nrow(sigres),"eQTL pairs with q<1e-5 and log2FC >1/<-1"))
print(paste(length(sigvariants), "unique subset of significant bolting-related variants are linked to"))
print(paste(length(sig_genes), "unique affected bolting regulators"))
save(list = c('output_eqtl', 'maybelinked', 'res', "sigres", "sig_genes", "sigvariants", "tophits"), file = "results/output_me.Rda")