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

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

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
  dplyr::select(sampleId, genomeId, chr, region_start = start, region_end = end, prophage_id) %>%
  # head(20) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    hgs = list(
      region_homology_groups(
        pandb = panOrgDb, genome = genomeId, chr = chr,
        start = region_start, end = region_end,
        cols = c("GID", "mRNA_id", "mRNA_key")
      )
    )
  )

phage_genes <- dplyr::select(proHgs, prophage_id, sampleId, hgs) %>%
  tidyr::unnest(cols = c(hgs)) %>%
  dplyr::rename(hgId = GID)

# write prophage region homology groups to a file
dplyr::group_by(phage_genes, prophage_id, sampleId) %>%
  dplyr::summarise(
    nHgs = length(hgId),
    hgs = paste(hgId, collapse = ";"),
    .groups = "drop"
  ) %>%
  dplyr::rename(id = prophage_id) %>%
  dplyr::ungroup() %>%
  dplyr::filter(nHgs > 0) %>%
  readr::write_tsv(
    file = confs$analysis$prophages$preprocessing$files$raw_prophage_hg
  )

readr::write_tsv(
  phage_genes,
  file = confs$analysis$prophages$preprocessing$files$raw_prophage_genes
)
