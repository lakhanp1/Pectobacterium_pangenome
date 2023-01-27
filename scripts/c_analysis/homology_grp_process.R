#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))

## 1) process homology groups to transpose the table where
## rows are genomes and columns are homology group identifiers
## 2) plot heatmap of homology group across pangenome

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name

analysisName <- "homology_groups"
outDir <- confs$analysis$association$dir
outPrefix <- file.path(outDir, analysisName)

################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

hgs <- suppressMessages(
  readr::read_csv(
    file = confs$data$pangenomes[[pangenome]]$db$gene_classification$GC.100.0$files$groups
  )
) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "Genome_", "")) %>% 
  dplyr::rename_with(.fn = ~tolower(.x))

g2hg <- dplyr::select(hgs, hg = homology_group_id, !!!sampleInfo$Genome) %>% 
  tidyr::pivot_longer(
    cols = -hg, names_to = "Genome", values_to = "nGenes"
  ) %>% 
  tidyr::pivot_wider(
    id_cols = Genome, names_from = hg, values_from = nGenes
  )

readr::write_tsv(
  x = g2hg, file = confs$analysis$homology_groups$files$groups
)

hgMeta <- dplyr::select(hgs, homology_group_id, class)

readr::write_tsv(x = hgMeta, file = confs$analysis$homology_groups$files$groups_meta)

################################################################################





