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

## visualize homology groups PAV for given homology group sets

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

analysisName <- "vir_lineage_prophages"

outDir <- file.path(confs$analysis$prophages$dir)
outPrefix <- file.path(outDir, analysisName)

setIds <- c(
  "g_188.vir_2", "g_400.vir_2", "g_399.vir_2", "g_448.vir_2", "g_189.vir_1", 
  "g_162.vir_4", "g_93.vir_1", "g_438.vir_3", "g_393.vir_1"
)

# add any homology group set file created to this vector
# `setIds` are searched in these files to identify respective homology groups
# for plotting
hgSetFiles <- c(
  confs$analysis$prophages$files$prophage_hg
)

################################################################################

treeMethod <- "ani_upgma" # ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

orgDb <- org.Pectobacterium.spp.pan.eg.db

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>% dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

################################################################################
# prepare homology groups data for plotting
hgSets <- NULL

for (st in hgSetFiles) {
  df <- suppressMessages(
    readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
  ) %>%
    dplyr::mutate(
      hgs = stringr::str_split(hgs, ";"),
      setSource = st
    ) %>% 
    dplyr::select(tidyselect::all_of(c("id", "hgs", "setSource")))
  
  hgSets <- dplyr::bind_rows(hgSets, df)
}

if (any(duplicated(hgSets$id))) {
  stop(
    "duplicate homology group set identifiers found in the data:",
    names(which(table(hgSets$id) > 1))
  )
}

# read prophage HGs stored locally
hgL <- dplyr::left_join(
  x = tibble::tibble(id = setIds), y = hgSets, by = "id"
) %>% 
  purrr::transpose() %>% 
  purrr::set_names(nm = purrr::map(., "id"))

################################################################################

htList <- species_key_heatmap(
  genomes = rawTree$tip.label, speciesInfo = sampleInfo
)

htList@heatmap_param$width <- unit(10, "cm")

for (id in names(hgL)) {
  
  # prepare homology group PAV matrix from pan.db
  hgMat <- homology_groups_mat(
    pandb = orgDb, type = "pav", groups = hgL[[id]]$hgs
  )
  
  hgMat <- hgMat[rawTree$tip.label, ]
  colnames(hgMat) <- paste("hg.", colnames(hgMat), sep = "")
  
  
  ht <- homology_group_heatmap(
    mat = hgMat, phy = rawTree,
    # width = unit(10, "cm"),
    name = id, column_title = id
  )
  
  ht@column_dend_param$cluster <- FALSE
  ht@column_names_param$show <- FALSE
  
  htList <- htList + ht
  
}


pdf(file = paste(outPrefix, ".hgs.pdf", sep = ""), width = 20, height = 12)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = hgL[[1]]$id,
  row_dend_side = "left",
  merge_legends = TRUE
)
dev.off()

################################################################################
