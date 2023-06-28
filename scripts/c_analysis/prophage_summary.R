#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))


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
panProphages <- dplyr::left_join(
  dplyr::select(sampleInfo, sampleId, Genome, SpeciesName, geo_loc_country), 
  y = prophageDf, by = "sampleId")

table(panProphages$checkv_quality, panProphages$SpeciesName)

dplyr::group_by(panProphages, sampleId, SpeciesName) %>% 
  dplyr::summarise(n = n_distinct(prophage_id, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(SpeciesName) %>% 
  skimr::skim()

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>% 
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>% 
  dplyr::arrange(desc(nHgs)) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName), by = "sampleId"
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(prophageDf, prophage_id, prophage_length,
                      completeness, checkv_quality, miuvig_quality),
    by = "prophage_id") %>% 
  dplyr::filter(prophage_length >= 5000)

################################################################################
# species wise prophage summary
proHgStats <- dplyr::group_by(proHgs, sampleId) %>% 
  dplyr::summarise(
    n_phages = n(),
    n_hgs = sum(nHgs),
    hgs = list(unlist(hgs)),
    .groups = "drop"
  ) %>% 
  dplyr::right_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName), by = "sampleId"
  ) %>% 
  tidyr::replace_na(replace = list(n_phages = 0, n_hgs = 0))

panVirSummary <- dplyr::summarise(
  proHgStats,
  n_genomes = n(),
  n_vir_sp = sum(n_phages),
  max_vir_per_g = max(n_phages),
  min_vir_per_g = min(n_phages),
  mean_vir_per_g = mean(n_phages),
  mean_vir_hgs = mean(n_hgs),
  median_vir_hgs = median(n_hgs),
  total_vir_hgs = sum(n_hgs),
  unique_vir_hgs = length(unique(unlist(hgs)))
) %>% 
  dplyr::mutate(SpeciesName = "Pangenome", .before = n_genomes)

spVirSummary <- dplyr::group_by(proHgStats, SpeciesName) %>% 
  dplyr::summarise(
    n_genomes = n(),
    n_vir_sp = sum(n_phages),
    max_vir_per_g = max(n_phages),
    min_vir_per_g = min(n_phages),
    mean_vir_per_g = mean(n_phages),
    mean_vir_hgs = mean(n_hgs),
    median_vir_hgs = median(n_hgs),
    total_vir_hgs = sum(n_hgs),
    unique_vir_hgs = length(unique(unlist(hgs)))
  ) %>%
  dplyr::bind_rows(panVirSummary) %>% 
  dplyr::arrange(desc(n_genomes)) %>% 
  dplyr::left_join(y = hgSummary, by = "SpeciesName")

readr::write_tsv(x = spVirSummary, file = confs$analysis$prophages$files$prophage_stats_species)

################################################################################

perGenomePhageInfo <- dplyr::left_join(
  proHgs, panProphages, by = c("contig_id", "prophage_id", "sampleId", "SpeciesName")
) %>% 
  dplyr::group_by(Genome, SpeciesName) %>% 
  dplyr::summarise(
    nProphages = n(),
    longest_phage = max(length),
    smallest_phage = min(length),
    mean_phage_len = mean(length),
    longest_nHg = max(nHgs),
    smallest_nHg = min(nHgs),
    .groups = "drop"
  ) %>% 
  dplyr::arrange(desc(nProphages))

readr::write_tsv(
  x = perGenomePhageInfo,
  file = confs$analysis$prophages$files$prophage_stats_genome
)



