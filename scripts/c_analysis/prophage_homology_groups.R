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

outDir <- confs$analysis$prophages$dir

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

prophageDf <- suppressMessages(
  readr::read_tsv(confs$data$prophages$files$data)
) %>%
  dplyr::rename(prophage_length = length) %>%
  dplyr::filter(viral_genes != 0) %>%
  dplyr::select(-SpeciesName, -Genome)

hgSummary <- suppressMessages(
  readr::read_tsv(confs$analysis$homology_groups$files$spp_group_stats)
) %>%
  tidyr::pivot_wider(
    id_cols = c(SpeciesName),
    names_from = class,
    values_from = count
  ) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    total_hgs = sum(core, accessory, unique, na.rm = T)
  )
################################################################################
panProphages <- dplyr::left_join(sampleInfo, prophageDf, by = "sampleId")

table(panProphages$checkv_quality, panProphages$SpeciesName)

dplyr::group_by(panProphages, sampleId, SpeciesName) %>%
  dplyr::summarise(n = n_distinct(prophage_id, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(SpeciesName) %>%
  skimr::skim()

# get homology groups for each prophage region
proHgs <- dplyr::filter(panProphages, !is.na(contig_id)) %>%
  # dplyr::select(sampleId, Genome, chr, start, end, contig_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    hgs = list(
      region_homology_groups(
        pandb = orgDb, genome = Genome, chr = chr, start = start, end = end
      )
    )
  )

# filter to keep high quality prophages
proHgs %>%
  dplyr::mutate(
    nHgs = length(hgs),
    hgs = paste(hgs, collapse = ";")
  ) %>%
  dplyr::filter(nHgs > 0) %>%
  dplyr::select(contig_id, prophage_id, sampleId, nHgs, hgs) %>%
  readr::write_tsv(
    file = confs$analysis$prophages$files$prophage_hg
  )

################################################################################
