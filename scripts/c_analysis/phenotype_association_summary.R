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

## plot PAV for the homology groups specific for a particular phenotype

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

phenotype <- "assay_FN"
treeMethod <- "ani_upgma"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name

outDir <- file.path(confs$analysis$association$dir, phenotype)
outPrefix <- file.path(outDir, phenotype)

################################################################################

if(!dir.exists(outDir)){
  dir.create(outDir)
}

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

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

## import homology group data
g2hg <- suppressMessages(
  readr::read_tsv(confs$analysis$homology_groups$files$groups)
) %>% 
  dplyr::mutate(Genome = as.character(Genome))

hgMeta <- suppressMessages(
  readr::read_tsv(confs$analysis$homology_groups$files$groups_meta)
)

phenoMeta <- suppressMessages(
  readr::read_csv(file = confs$analysis$association$files$phenotypes)
) %>% 
  dplyr::select(Genome, all_of(phenotype)) %>% 
  dplyr::filter(!is.na(!!sym(phenotype))) %>% 
  dplyr::group_by(!!sym(phenotype))

## process phenotype association results
res <- suppressMessages(
  readr::read_csv(
    file = confs$data$pangenomes[[pangenome]]$db$gene_classification$phenotypes[[phenotype]]$files$association,
    comment = "#"
  )
) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
  dplyr::rename_with(.fn = ~tolower(.x)) %>% 
  dplyr::filter(phenotype == phenotype, fisher_exact_p_value != "No test")


################################################################################
## homology groups specific for a phenotype
phenoSpecific <- dplyr::filter(
  res, 
  phenotype_members_absent == 0,
  other_phenotype_members_present == 0
)

dplyr::select(phenoSpecific, phenotype, homology_group_id) %>% 
  dplyr::group_by(phenotype) %>% 
  dplyr::summarize(
    homology_group_id = stringr::str_c(homology_group_id, collapse = ",")
  ) %>% 
  readr::write_tsv(
    file = confs$analysis$association$files$pheno_specific_groups,
    append = TRUE
  )

## remove duplicates and rewrite the data
suppressMessages(
  readr::read_tsv(
    confs$analysis$association$files$pheno_specific_groups, col_names = FALSE,
    col_types = "cc"
  )
) %>% 
  dplyr::distinct() %>% 
  readr::write_tsv(
    file = confs$analysis$association$files$pheno_specific_groups,
    col_names = FALSE
  )

################################################################################
#' Plot homology group heatmap with provided clustering
#'
#' @param mat homology group matrix
#' @param phy A phylo object
#' @param metadata metadata for genomes. Must have SpeciesName and Genome column
#' @param width Vector of length 2 for heatmap widths
#'
#' @return ComplexHeatmap with layout: species key heatmap | ANI heatmap
#' @export
#'
#' @examples
homology_group_heatmap <- function(mat, phy, metadata, width){
  
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    # https://github.com/jokergoo/ComplexHeatmap/issues/949
    all(rownames(mat) == phy$tip.label),
    all(rownames(mat) %in% metadata$Genome),
    any(class(phy) == "phylo"),
    length(width) == 2,
    all(is.numeric(width))
  )
  
  ## homology group heatmap
  ht_hg <- ComplexHeatmap::Heatmap(
    matrix = hgMat,
    name = "hg",
    col = circlize::colorRamp2(
      breaks = seq(0, min(3, max(hgMat))),
      colors = viridisLite::viridis(n = min(3, max(hgMat))+1, option = "B")
    ),
    # heatmap_legend_param = list(legend_height = unit(3, "cm")),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE, show_column_names = FALSE,
    column_title = "Homology groups",
    width = unit(width[2], "cm")
  )
  
  
  # draw(ht_hg)
  
  ## species name key heatmap
  speciesMat <- tibble::tibble(
    Genome = rownames(mat)
  ) %>% 
    dplyr::left_join(y = dplyr::select(metadata, Genome, SpeciesName), by = "Genome") %>% 
    dplyr::mutate(sp = 1) %>% 
    tidyr::pivot_wider(
      id_cols = Genome, names_from = SpeciesName,
      values_from = sp, values_fill = 0, names_sort = TRUE
    ) %>% 
    tibble::column_to_rownames(var = "Genome") %>% 
    as.matrix()
  
  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == phy$tip.label))
  
  ht_species <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = c("1" = "black", "0" = "white"),
    # cluster_rows = as.hclust.phylo(phy), row_dend_reorder = FALSE,
    cluster_columns = FALSE,
    # column_order = speciesOrder,
    column_split = 1:ncol(speciesMat), cluster_column_slices = FALSE,
    border = TRUE, column_gap = unit(0, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_title = "Species key",
    width = unit(width[1], "cm")
  )
  
  # draw(ht_species)
  
  htList <- ht_species + ht_hg
  
  return(htList)
  
}

################################################################################
hgMat <- dplyr::select(g2hg, Genome, as.character(phenoSpecific$homology_group_id)) %>% 
  # dplyr::left_join(y = phenoMeta, by = "Genome") %>% 
  # dplyr::relocate(!!phenotype, .after = Genome)
  tibble::column_to_rownames(var = "Genome") %>% 
  as.matrix()

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo, width = c(10,12)
)

png(filename = paste(outPrefix, ".pheno_hg_association.png", sep = ""), width = 5000, height = 2800, res = 350)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
  
)
dev.off()


# nodeOfInterest <- dplyr::filter(sampleInfo, SpeciesName == "P. brasiliense") %>% 
#   dplyr::pull(Genome)
# 
# clade <- ape::getMRCA(phy = rawTree, tip = nodeOfInterest)
# subTree <- ape::extract.clade(phy = rawTree, node = clade)

