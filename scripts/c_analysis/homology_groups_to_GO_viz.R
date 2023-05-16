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

## visualize homology groups for GO terms

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

analysisName <- "virulence_GO"

goIds <- c(
  "GO:0046718", "GO:0044826", "GO:0019069", "GO:0044659", "GO:0075713",
  "GO:0044660", "GO:0019076", "GO:0019058"
)

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

hgs <- AnnotationDbi::select(
  x = orgDb, keys = goIds, keytype = "GOALL", columns = "GID"
) %>% 
  dplyr::mutate(count = 1) %>% 
  dplyr::rename(hg_id = GID) %>% 
  tidyr::pivot_wider(
    id_cols = hg_id, names_from = GOALL, values_from = count, values_fill = 0
  )


################################################################################
# prepare homology group PAV matrix from pan.db
hgMat <- homology_groups_mat(pandb = orgDb, type = "cnv", groups = hgs$hg_id)

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo,
  width = c(9, 20)
)

# GO term top annotation
goMat <- tibble::column_to_rownames(hgs, var = "hg_id") %>% as.matrix()
go_lables <- AnnotationDbi::mapIds(
  x = GO.db, keys = colnames(goMat), keytype = "GOID", column = "TERM"
)

colnames(goMat) <- go_lables

an_go <- ComplexHeatmap::HeatmapAnnotation(
  go = goMat,
  name = "go",
  which = "column",
  col = list(go = c("1" = "black", "0" = "white")),
  annotation_name_side = "left",
  annotation_name_gp = gpar(fontsize = 12),
  height = unit(5, "cm")
)

# add top annotation to the main heatmap
htList@ht_list$hg@top_annotation <- an_go

htList@ht_list$hg@column_names_param$show <- FALSE
htList@ht_list$hg@column_dend_param$cluster <- FALSE



pdf(file = paste(outPrefix, ".homology_grps.pdf", sep = ""), width = 14, height = 10)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
  merge_legends = TRUE
)
dev.off()

################################################################################


