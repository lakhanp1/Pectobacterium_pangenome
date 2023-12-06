#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# extract homology groups for plasmid regions

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

plasmidDf <- suppressMessages(
  readr::read_tsv(confs$data$plasmids$files$data)
) %>%
  dplyr::select(-SpeciesName, -genomeId, plasmid_length = length)

################################################################################
panPlasmids <- dplyr::left_join(sampleInfo, plasmidDf, by = "sampleId")

dplyr::group_by(panPlasmids, sampleId, SpeciesName) %>%
  dplyr::summarise(n = n_distinct(plasmid_id, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(SpeciesName) %>%
  skimr::skim()

# get homology groups for each plasmid region
plasmidHgs <- dplyr::filter(panPlasmids, !is.na(plasmid_id)) %>% 
  dplyr::select(sampleId, genomeId, chr, plasmid_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    hgs = list(
      region_homology_groups(
        pandb = orgDb, genome = genomeId, chr = chr
      )
    )
  )

# filter to remove prophages without any hgs 
plasmidHgs %>%
  dplyr::mutate(
    nHgs = length(hgs),
    hgs = paste(hgs, collapse = ";")
  ) %>%
  dplyr::filter(nHgs > 0) %>%
  dplyr::select(plasmid_id, sampleId, nHgs, hgs) %>%
  readr::write_tsv(
    file = confs$analysis$plasmids$preprocessing$files$raw_plasmid_hg
  )

################################################################################
