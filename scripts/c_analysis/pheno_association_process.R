#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))

## save phenotype specific groups in a file
##**********
##*replace this script later by a function
##**********

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/association_analysis.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
pangenomeConf <- confs$data$pangenomes[[pangenome]]

pangenomeConf$db$gene_classification$phenotypes$dir

################################################################################

phnConf <- suppressMessages(
  readr::read_tsv(
    file = pangenomeConf$analysis_confs$files$phenotype_association, 
    col_types = "ccccc"
  )
)

phn <- "assay_FN"

## extract and save phenotype specific groups
for (phn in phnConf$name) {
  cat("phenotype", phn, "\n")
  phenotype_specific_groups(
    phenotype = phn, panConf = pangenomeConf,
   save = confs$analysis$association$files$pheno_specific_groups 
  )
}

