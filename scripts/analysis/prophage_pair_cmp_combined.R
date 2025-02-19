#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(configr))

# combine prophage comparison metrics (ANI, MASH, syntenic Jaccard) into one table

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

outDir <- confs$analysis$prophages$preprocessing$path

panOrgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################
# prophage clusters
regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
) %>%
  dplyr::select(
    id = prophage_id, grp = phage_grp, genomeId,
    host = SpeciesName, loc = geo_loc_country
  )

# MASH
mashDist <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$mash_dist)
) %>%
  tidyr::pivot_longer(!prophage_id, names_to = "pro2", values_to = "mash") %>%
  dplyr::rename(pro1 = prophage_id)

# ANI
aniDf <- suppressMessages(readr::read_tsv(
  file = confs$data$prophages$files$ani,
  col_names = c("pro1", "pro2", "ani", "mapped", "total")
)) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(pro1, pro2),
      .fns = ~ stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1")
    ),
    dist = 1 - (ani / 100)
  ) %>%
  dplyr::select(pro1, pro2, ani)

# HG syntenic Jaccard index
simDf <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$pair_comparison)
) %>%
  dplyr::select(phage1, phage2, jaccardIndex, nSharedHgs, nSyntenicHgs)

phageCmpDf <- dplyr::bind_rows(
  dplyr::rename(simDf, pro1 = phage1, pro2 = phage2),
  dplyr::rename(simDf, pro2 = phage1, pro1 = phage2)
)

################################################################################
# combine data
pairwiseDf <- tidyr::expand_grid(pro1 = regionClusters$id, pro2 = regionClusters$id) %>%
  dplyr::left_join(
    y = dplyr::rename_with(regionClusters, .fn = ~ paste("pro1_", .x, sep = "")),
    by = c("pro1" = "pro1_id")
  ) %>%
  dplyr::left_join(
    y = dplyr::rename_with(regionClusters, .fn = ~ paste("pro2_", .x, sep = "")),
    by = c("pro2" = "pro2_id")
  ) %>%
  dplyr::filter(pro1_grp == pro2_grp, pro1 != pro2)

phageCmpDf <- dplyr::left_join(
  x = pairwiseDf, y = mashDist, by = c("pro1", "pro2")
) %>%
  dplyr::left_join(y = aniDf, by = c("pro1", "pro2")) %>%
  dplyr::left_join(y = phageCmpDf, by = c("pro1", "pro2"))

readr::write_tsv(
  phageCmpDf,
  file = paste(confs$analysis$prophages$path, "/intercluster_pairwise_cmp.tab", sep = "")
)
