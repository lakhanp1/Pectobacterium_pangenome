#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# extract homology groups for prophage regions

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

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId 
)

prophageDf <- suppressMessages(
  readr::read_tsv(confs$data$prophages$files$data)
) %>%
  dplyr::filter(viral_genes != 0) %>%
  dplyr::select(-SpeciesName, -genomeId)

################################################################################
panProphages <- dplyr::left_join(sampleInfo, prophageDf, by = "sampleId")

table(panProphages$checkv_quality, panProphages$SpeciesName)

dplyr::group_by(panProphages, sampleId, SpeciesName) %>%
  dplyr::summarise(n = n_distinct(prophage_id, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(SpeciesName) %>%
  skimr::skim()

# get homology groups for each prophage region
proHgs <- dplyr::filter(panProphages, !is.na(contig_id)) %>%
  dplyr::select(sampleId, genomeId, chr, start, end, prophage_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    hgs = list(
      region_homology_groups(
        pandb = orgDb, genome = genomeId, chr = chr, start = start, end = end
      )
    )
  )

# filter to remove prophages without any hgs 
proHgs %>%
  dplyr::mutate(
    nHgs = length(hgs),
    hgs = paste(hgs, collapse = ";")
  ) %>%
  dplyr::filter(nHgs > 0) %>%
  dplyr::select(id = prophage_id, sampleId, nHgs, hgs) %>%
  readr::write_tsv(
    file = confs$analysis$prophages$preprocessing$files$raw_prophage_hg
  )

################################################################################
