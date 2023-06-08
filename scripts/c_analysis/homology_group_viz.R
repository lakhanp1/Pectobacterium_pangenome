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
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

## visualize homology groups PAV 

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/association_analysis.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

analysisName <- "prophage_hgs"

hgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>%
  dplyr::filter(!is.na(hgs)) %>%
  dplyr::filter(
    sampleId %in% c("JGAR-100719-1", "NAK238", "NAK716", "NAK700", "GCF_000754695.1_ASM75469v1")
    # sampleId == "JGAR-100719-1"
  ) %>%
  # dplyr::pull(hgs) %>%
  stringr::str_split(";") %>%
  unlist() %>% unique()

treeMethod <- "ani_upgma"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- file.path(confs$analysis$homology_groups$dir)
outPrefix <- file.path(outDir, analysisName)

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)


################################################################################
# prepare homology group PAV matrix from pan.db
hgMat <- homology_groups_mat(pandb = orgDb, type = "pav", groups = hgs)

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo,
  width = c(10, 20)
)

htList@ht_list$hg@column_dend_param$cluster <- FALSE 
htList@ht_list$hg@column_names_param$show <- FALSE

pdf(file = paste(outPrefix, ".homology_grps.pdf", sep = ""), width = 20, height = 9)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
  merge_legends = TRUE
)
dev.off()

################################################################################


