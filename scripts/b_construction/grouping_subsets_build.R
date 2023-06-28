suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ape))

## perform optimal grouping on subset of the data:
## 1) 20 clades of UPGMA tree
## 2) only type strains

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- here::here("analysis", "03_pangenome_pecto_v2", "subset_optimal_group")

################################################################################

if (!dir.exists(outDir)) {
  dir.create(outDir, recursive = T)
}

metadata <- suppressMessages(readr::read_tsv(confs$data$reference_data$files$metadata)) %>%
  dplyr::select(sampleId, length, GC_per, n_contigs, N50, L50)

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) %>%
  dplyr::left_join(y = metadata, by = c("sampleId"))

treeUpgma <- ape::read.tree(file = confs$analysis$phylogeny$ani_upgma$files$tree)
################################################################################
## type strains, with best N50 in case of multiple type strains per species
typeStrains <- dplyr::filter(sampleInfo, !is.na(type_material)) %>%
  dplyr::group_by(SpeciesName) %>%
  dplyr::arrange(desc(N50)) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

################################################################################
# cut tree into 20 clusters and (select 1 random representative from each cluster) * 100 times
treeSets <- tibble::enframe(
  x = cutree(tree = as.hclust.phylo(treeUpgma), k = 20),
  name = "id", value = "group"
) %>%
  dplyr::add_count(group, name = "count") %>%
  dplyr::arrange(count) %>%
  dplyr::group_by(group) %>%
  dplyr::slice_sample(n = 100, replace = TRUE) %>%
  dplyr::mutate(
    set = 1:n(),
    setId = paste("rand", sprintf(fmt = "%03d", set), sep = "_")
  ) %>%
  dplyr::ungroup() %>%
  dplyr::left_join(y = dplyr::select(sampleInfo, sampleId, Genome), by = c("id" = "Genome"))


optGrpSets <- dplyr::group_by(treeSets, setId) %>%
  dplyr::summarise(genomes = paste(id, collapse = ",")) %>%
  dplyr::bind_rows(
    tibble(
      setId = "typeStrain", genomes = paste(typeStrains$Genome, collapse = ",")
    )
  ) %>%
  dplyr::mutate(
    setId = forcats::fct_relevel(setId, "typeStrain")
  ) %>%
  dplyr::arrange(setId)

readr::write_tsv(
  x = optGrpSets, col_names = F,
  file = file.path(outDir, "subsets_conf.tab")
)
