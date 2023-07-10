#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(visNetwork))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggraph))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

# visualize representative prophages in pangenome

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

treeMethod <- "ani_upgma"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- confs$analysis$prophages$dir
outPrefix <- ""

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>% 
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName, Genome, N50), by = "sampleId"
  ) 

prophageDf <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>% 
  dplyr::select(prophage_id, prophage_length = length, completeness)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

phageRelations <- suppressMessages(
  readr::read_tsv(file = confs$analysis$prophages$files$network)
) %>% 
  dplyr::left_join(y = prophageDf, by = c("prophage_id"))

phageTree <- ape::read.tree(confs$analysis$prophages$files$mash_upgma)
################################################################################
# prophage homology groups on pangenome phylogeny
proReps <- dplyr::filter(phageRelations, nodeType != "child") %>% 
  dplyr::select(prophage_id, nodeType) %>% 
  dplyr::left_join(proHgs, by = "prophage_id") %>% 
  dplyr::select(-contig_id)

length(unique(unlist(proReps$hgs)))

# function to get homology group counts in each genome of a pangenome
homology_groups_genome_counts <- function(h, pandb){
  # ****
  # replace this function by searching for a tandem match of homology groups
  # ****
  suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = unique(h),
      columns = c("genome")
    )
  ) %>% 
    dplyr::group_by(genome) %>% 
    dplyr::summarise(n = n())
}

proHgPerGenome <- dplyr::select(proReps, prophage_id, hgs, nHgs) %>% 
  # dplyr::slice_head(n = 10) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(
    genomeCounts = list(homology_groups_genome_counts(h = hgs, pandb = orgDb)
    )
  ) %>% 
  dplyr::select(-hgs) %>% 
  tidyr::unnest(cols = genomeCounts) %>% 
  dplyr::mutate(fraction = n/nHgs)

proHgCounts <- tidyr::pivot_wider(
  data = proHgPerGenome,
  id_cols = c(genome),
  names_from = prophage_id,
  values_from = fraction,
  values_fill = 0,
) %>% 
  dplyr::select(genome, all_of(phageTree$tip.label))


countMat <- tibble::tibble(genome = rawTree$tip.label) %>% 
  dplyr::left_join(y = proHgCounts, by = c("genome")) %>% 
  tibble::column_to_rownames(var = "genome") %>% 
  as.matrix()

hist(matrixStats::colMaxs(countMat))
hist(countMat)

# species key heatmap
ht_species <- species_key_heatmap(
  genomes = rawTree$tip.label, speciesInfo = sampleInfo
)

prophageAnDf <- tibble::tibble(prophage_id = colnames(countMat)) %>% 
  dplyr::left_join(y = phageRelations, by = "prophage_id")

## ensure the row order is same: this is because of a bug in ComplexHeatmap
stopifnot(
  all(rownames(countMat) == rawTree$tip.label),
  all(rownames(speciesMat) == rawTree$tip.label)
)

# bottom annotation: number of HGs in prophage, prophage completeness
an_bar <- ComplexHeatmap::HeatmapAnnotation(
  nHgs = ComplexHeatmap::anno_barplot(
    x = prophageAnDf$nHgs,
    axis_param = list(side = "right")
  ),
  completeness = ComplexHeatmap::anno_barplot(
    x = prophageAnDf$completeness,
    gp = gpar(col = "blue"),
    axis_param = list(side = "right")
  ),
  which = "column",
  annotation_name_side = "right",
  annotation_name_rot = 0,
  annotation_name_gp = gpar(fontsize = 12),
  annotation_height = unit(c(2, 0.5), "cm")
)

# prophage homology group count heatmap
ht_count <- ComplexHeatmap::Heatmap(
  matrix = countMat,
  name = "hg_count",
  col = viridisLite::viridis(n = 11, option = "B"),
  bottom_annotation = an_bar,
  cluster_rows = ape::as.hclust.phylo(rawTree), row_dend_reorder = FALSE,
  row_dend_width = unit(3, "cm"), show_row_dend = TRUE,
  cluster_columns = ape::as.hclust.phylo(phageTree), column_dend_reorder = FALSE,
  column_dend_height = unit(3, "cm"), show_column_dend = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  column_title = "Prophage homology groups / genome",
  width = unit(20, "cm")
)

htList <- ht_species + ht_count

pdf(file = paste(outDir, "/prophage_pangenome_hgs.pdf", sep = ""), width = 15, height = 9)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg_count",
  row_dend_side = "left",
  merge_legends = TRUE
)
dev.off()

