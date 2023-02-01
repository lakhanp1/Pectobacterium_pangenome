#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))

## prepare configuration files for running gene_classification with phenotype
## generated for various clades

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name

################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

cladeGrps <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$files$clade_copare)
) %>% 
  purrr::transpose() %>% 
  purrr::set_names(nm = purrr::map(., "name"))

# ## test
# cmp <- cladeGrps$NL_clade_vs_rest_Pbr
# treeName <- cmp$tree
# clade_comparison_confs(
#   file = confs$analysis$phylogeny[[treeName]]$files$tree_rooted,
#   node = cmp$node,
#   ancestorClade = cmp$clade,
#   name = cmp$name,
#   category = cmp$name
# )
# ######################

## get the genome IDs inside a clade being compared, its phenotype df
## and optionally genome ids for clade against the comparison will be made
cladeCmpList <- purrr::map(
  .x = cladeGrps,
  .f = function(cmp){
    treeName <- cmp$tree
    
    clade_comparison_confs(
      file = confs$analysis$phylogeny[[treeName]]$files$tree_rooted,
      node = cmp$node,
      ancestorClade = cmp$clade,
      name = cmp$name,
      category = cmp$name
    )
    
  }
)

################################################################################
## write updated phenotype file with clade of interest as phenotypes
cladePhenotypes <- dplyr::select(sampleInfo, Genome, everything(), -nodeLabs)

for (cmp in cladeCmpList) {
  # print(cmp$name)
  cladePhenotypes <- dplyr::left_join(
    cladePhenotypes, cmp$pheno, by = "Genome"
  )
}

readr::write_csv(
  cladePhenotypes,
  file = confs$analysis$phylogeny$files$clade_phenotypes
)


## write `pantools gene_classification` configuration
## format: <name>\t<--phenotype="xyz" [--include="1,3,5,15,..N"]>\t<compare_genomes>\t<against_genomes>\t
purrr::map_dfr(
  cladeCmpList,
  .f = function(x){
    return(
      list(
        name = x$name, include = x$includeSet,
        compare = x$compare, against = x$against
      )
    )
  }
) %>% 
  dplyr::mutate(
    phenotypeArg = dplyr::if_else(
      condition = is.na(include) | include == "",
      true = stringr::str_c("--phenotype=", name, sep = ""), 
      false = stringr::str_c("--phenotype=", name, " --include=", include, sep = "")
    )
  ) %>% 
  dplyr::select(name, phenotypeArg, compare, against, include) %>% 
  readr::write_tsv(
    file = confs$data$analysis_confs$files$phenotype_association,
    col_names = FALSE, na = ""
  )


################################################################################


