#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))

## generate input files for pyseer association analysis of virulence in Pbr isolates

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

analysisName <- "pyseer_clades"
outDir <- file.path(confs$analysis$association$dir, analysisName)
outPrefix <- file.path(outDir, analysisName)
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

if (!dir.exists(outDir)) {
  dir.create(outDir)
}
################################################################################
## clade as phenotypes
pheno <- suppressMessages(
  readr::read_csv(panConf$analysis_confs$files$clade_phenotypes)
) %>%
  dplyr::mutate(
    Genome = as.character(Genome),
    dplyr::across(
      .cols = -Genome,
      .fns = ~ dplyr::if_else(.x == "N", true = 0, false = 1, missing = NULL)
    )
  )

pheno %<>% dplyr::left_join(
  y = dplyr::select(sampleInfo, samples = sampleId, Genome, SpeciesName), by = "Genome"
) %>%
  dplyr::select(
    samples, Genome, SpeciesName, everything()
  )

readr::write_tsv(
  pheno,
  file = paste(outDir, "/clade.phenotypes.tab", sep = "")
)

## homology groups
hgs <- homology_groups_extract(
  file = panConf$db$gene_classification$GC.100.0$files$groups,
  genomes = sampleInfo$Genome,
  groups <- "accessory", pav = TRUE
)

dplyr::select(hgs, -class, Gene = hg, !!!purrr::map_chr(sampleInfoList, "Genome")) %>%
  dplyr::mutate(Gene = paste("hg_", Gene, sep = "")) %>%
  readr::write_tsv(
    file = paste(outDir, "/clade.accessory_PAV.tab", sep = "")
  )
