## jbrowse.R
# Miguel Alburo
# 03/05/24
# Checking whether variants are nonsynonymous.

# Packages
library(tidyverse)
# library(JBrowseR)
#
# # Create a new JBrowse instance
# JBrowseR("ViewHg19", location = "10:29,838,737..29,838,819")

source("https://bioconductor.org/biocLite.R")
biocLite("TxDb.Athaliana.BioMart.plantsmart28")
library(TxDb.Athaliana.BioMart.plantsmart28)