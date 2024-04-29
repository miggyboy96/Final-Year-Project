##
# Miguel Alburo
# 12/03/2024
#

## Load packages
library(tidyverse)
library(readxl)
source(file = 'code/functions.R')
library(BioNERO)

# Load data
grn_data <- readxl::read_excel(path = 'data/raw/media-10.xlsx')
grn_data <- as.data.frame(grn_data) %>%
  select("regulatory gene", "target gene", weight) #%>%
  #slice(which(! grn_data$`target gene` %in% grn_data$`regulatory gene`))

## Print GRN
a <- plot_grn(
  edgelist_grn = grn_data,
  show_labels = "tophubs",
  top_n_hubs = 44,
  layout = igraph::with_kk,
  arrow.gap = 0.01,
  ranked = TRUE,
  curvature = 0.1,
  interactive = FALSE,
  dim_interactive = c(600, 600)
)