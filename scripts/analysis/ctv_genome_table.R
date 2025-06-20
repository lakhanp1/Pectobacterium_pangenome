suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))


rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- paste(confs$analysis$ctv$data$path)

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) %>%
  dplyr::select(genomeId, SpeciesName)

ctvCluster <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
) %>%
  dplyr::filter(phage_grp == "phage_grp_1") %>%
  dplyr::select(prophage_id, fragments, nFragments, nHg, genomeId, completeness)

rawProphages <- suppressMessages(
    readr::read_tsv(confs$data$prophages$files$data)
  )

prophagePos <- dplyr::select(ctvCluster, genomeId, prophage_id, fragments) %>%
  dplyr::mutate(
    fragments = stringr::str_split(fragments, pattern = ";")
  ) %>%
  tidyr::unnest(fragments) %>%
  dplyr::left_join(
    y = dplyr::select(rawProphages, prophage_id, chr, start, end), by = "prophage_id") %>%
  dplyr::mutate(
    pos = paste(chr, ":", start, "-", end, sep = "")
  ) %>%
  dplyr::group_by(genomeId, prophage_id) %>%
  dplyr::summarise(
    fragments = paste(fragments, collapse = ";"),
    pos = paste(pos, collapse = ";"),
    .groups = "drop"
  )


ctvCluster <- dplyr::select(ctvCluster, -fragments) %>%
  dplyr::left_join(prophagePos, by = c("genomeId", "prophage_id")) %>%
  dplyr::relocate(fragments, .after = prophage_id)

sampleInfo <- dplyr::left_join(
  x = sampleInfo, y = ctvCluster, by = "genomeId"
)

readr::write_tsv(
  sampleInfo,
  file = file.path(outDir, "ctv_genome_table.tsv")
)
