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

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
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

# phn <- "assay_FN"
for (phn in phnConf$name) {
  cat("phenotype", phn, "\n")
  
  file_phenoAsso <- file.path(
    pangenomeConf$db$gene_classification$phenotypes$dir, phn, "phenotype_association.csv"
  )
  
  ## process phenotype association results
  res <- suppressMessages(
    readr::read_csv(file = file_phenoAsso, comment = "#")
  ) %>% 
    dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
    dplyr::rename_with(.fn = ~tolower(.x)) %>% 
    dplyr::filter(phenotype == phn, fisher_exact_p_value != "No test")
  
  ## homology groups specific for a phenotype
  phenoSpecific <- dplyr::filter(
    res, 
    phenotype_members_absent == 0,
    other_phenotype_members_present == 0
  )
  
  if(nrow(phenoSpecific) > 0){
    dplyr::select(phenoSpecific, phenotype, homology_group_id) %>% 
      dplyr::group_by(phenotype) %>% 
      dplyr::summarize(
        homology_group_id = stringr::str_c(homology_group_id, collapse = ",")
      ) %>% 
      readr::write_tsv(
        file = confs$analysis$association$files$pheno_specific_groups,
        append = TRUE
      )
  }
}


## remove duplicates and rewrite the data
suppressMessages(
  readr::read_tsv(
    confs$analysis$association$files$pheno_specific_groups, col_names = FALSE,
    col_types = "cc"
  )
) %>% 
  dplyr::distinct() %>% 
  readr::write_tsv(
    file = confs$analysis$association$files$pheno_specific_groups,
    col_names = FALSE
  )



