#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))

## prepare configuration files for running gene_classification with phenotype
## generated for various clades

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
################################################################################

assemblyMeta <- suppressMessages(readr::read_tsv(confs$data$reference_data$files$metadata)) %>% 
  dplyr::select(sampleId, length,N50, N90, L50, L90)

sampleInfo <- get_metadata(file = panConf$files$metadata) %>% 
  dplyr::left_join(y = assemblyMeta, by = "sampleId")

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

cladeGrps <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$files$clade_copare, comment = "#")
) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = c(compare, background, bgExcludeNode, bgExcludeTip),
      .fns = ~stringr::str_split(string = .x, pattern = ",")
    )
  ) %>% 
  purrr::transpose() %>% 
  purrr::set_names(nm = purrr::map(., "name"))

# ## test
# cmp <- cladeGrps$versatile_clade2
# treeName <- cmp$tree
# 
# clade_comparison_confs(
#   file = confs$analysis$phylogeny[[treeName]]$files$tree_rooted,
#   node = cmp$compare,
#   type = cmp$compareType,
#   against = cmp$background,
#   name = cmp$name,
#   category = cmp$name,
#   excludeNode = cmp$bgExcludeNode,
#   excludeTips = cmp$bgExcludeTip
# )
# ######################

## get the genome IDs inside a clade being compared, its phenotype df
## and optionally genome ids for clade against the comparison will be made
cladeCmpList <- purrr::map(
  .x = cladeGrps,
  .f = function(cmp){
    treeName <- cmp$tree
    
    clade_comparison_confs(
      tree = confs$analysis$phylogeny[[treeName]]$files$tree_rooted,
      node = cmp$compare,
      type = cmp$compareType,
      against = cmp$background,
      name = cmp$name,
      category = cmp$name,
      excludeNode = cmp$bgExcludeNode,
      excludeTips = cmp$bgExcludeTip
    )
    
  }
)

################################################################################
## write updated phenotype file with clade of interest as phenotypes
cladePhenotypes <- dplyr::select(sampleInfo, Genome)

for (cmp in cladeCmpList) {
  # print(cmp$name)
  cladePhenotypes <- dplyr::left_join(
    cladePhenotypes, cmp$pheno, by = "Genome"
  )
}

readr::write_csv(
  cladePhenotypes,
  file = panConf$analysis_confs$files$clade_phenotypes
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
  dplyr::mutate(compare = stringr::str_split(compare, ",")) %>% 
  tidyr::unnest(cols = compare) %>% 
  dplyr::left_join(y = dplyr::select(sampleInfo, Genome, N50, length),
                   by = c("compare" = "Genome")) %>% 
  dplyr::group_by(name) %>% 
  dplyr::arrange(desc(N50), .by_group = TRUE) %>% 
  dplyr::mutate(compare = stringi::stri_flatten(compare, collapse = ",")) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>%
  dplyr::select(-N50, -length) %>% 
  dplyr::mutate(
    phenotypeArg = dplyr::if_else(
      condition = is.na(include) | include == "",
      true = stringr::str_c("--phenotype=", name, sep = ""), 
      false = stringr::str_c("--phenotype=", name, " --include=", include, sep = "")
    ),
    name = forcats::fct_relevel(
      name, purrr::map_chr(cladeGrps, "name") %>% unname
    )
  ) %>% 
  dplyr::arrange(name) %>% 
  dplyr::select(name, phenotypeArg, compare, against, include) %>% 
  readr::write_tsv(
    file = panConf$analysis_confs$files$clade_association,
    na = ""
  )


################################################################################
## save config for species wise analysis 
spGenomes <- dplyr::group_by(sampleInfo, SpeciesName) %>% 
  dplyr::summarise(
    count = n(),
    genomes = paste(Genome, collapse = ",")
  ) %>% 
  dplyr::arrange(desc(count)) %>% 
  dplyr::mutate(
    SpeciesName = stringr::str_replace_all(
      string = SpeciesName, pattern = "\\W+", replacement = "_"
    )
  )

readr::write_tsv(
  spGenomes, file = panConf$analysis_confs$files$species_genomes
)
