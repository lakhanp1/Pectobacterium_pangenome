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

## plot PAV for the homology groups specific for a particular phenotype

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/GO_enrichment/enrichment_functions.R")
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

phenotype <- "assay_FN"
treeMethod <- "ani_upgma"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- file.path(confs$analysis$association$dir, phenotype)
outPrefix <- file.path(outDir, phenotype)

## sequence info file for mRNAs across all genomes belonging to homology groups 
## specific to a phenotype of interest
file_associatedSeqInfo <- paste(outPrefix, ".pheno_specific.seq_info.txt", sep = "")
orgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

if(!dir.exists(outDir)){
  dir.create(outDir)
}

sampleInfo <- get_metadata(file = panConf$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## genomes that are used for phenotype association analysis
associatedGenomes <- get_phenotype_association_genomes(
  phenotype = phenotype,
  confFile = panConf$analysis_confs$files$clade_association
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

## get phenotype specific groups
specGrps <- phenotype_specific_groups(
  phenotype = phenotype, panConf = panConf
) %>% as.character()


## get the best genome
hgSeqInfo <- AnnotationDbi::select(
  x = orgDb,
  keys = specGrps,
  columns = c(
    "genome", "chr", "chr_id", "chr_name", "start", "end", "strand", "mRNA_id", "mRNA_key"
  )
) %>% 
  dplyr::rename(Genome = genome, homology_group_id = GID) %>% 
  dplyr::mutate(
    dplyr::across(.cols = c(chr, start, end), .fns = as.integer)
  ) %>% 
  tibble::as_tibble()


genomeSummary <- dplyr::add_count(hgSeqInfo, Genome, name = "nHg") %>% 
  dplyr::group_by(Genome, nHg) %>% 
  dplyr::summarise(nChr = length(unique(chr)), .groups = "drop") %>% 
  dplyr::left_join(
    dplyr::select(sampleInfo, Genome, SpeciesName, geo_loc_country, host, isolation_source,
                  env_broad_scale, collection_year),
    by = "Genome") %>% 
  dplyr::arrange(desc(nHg), nChr)

readr::write_tsv(
  x = genomeSummary,
  file = paste(outPrefix, ".specific_hgs.pangenome_summary.tab", sep = "")
)

bestGenome <- associatedGenomes$compare[1]

hgAnnotation <- dplyr::filter(hgSeqInfo, Genome == !!bestGenome) %>% 
  dplyr::arrange(chr, start) %>% 
  dplyr::mutate(
    label = paste(mRNA_id, "| ", chr, ":", start, sep = "")
  )

################################################################################
# functional annotation for homology groups
# COG
grpCog <- AnnotationDbi::select(
  x = orgDb,
  keys = hgAnnotation$homology_group_id,
  columns = c("GID", "COG_description", "COG_category")
) %>% 
  dplyr::rename(homology_group_id = GID) %>% 
  dplyr::left_join(y = hgAnnotation, by = "homology_group_id") %>% 
  dplyr::select(homology_group_id, Genome, chr, start, end, strand, mRNA_id,
                starts_with("COG"))

# topGO
topGoDf <- topGO_enrichment(genes = specGrps, orgdb = orgDb)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "COG")
openxlsx::writeData(
  wb = wb, sheet = 1, x = grpCog,
  startCol = 1, startRow = 1, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 2, firstActiveCol = 2)

openxlsx::addWorksheet(wb = wb, sheetName = "topGO")
openxlsx::writeData(
  wb = wb, sheet = 2, x = topGoDf,
  startCol = 1, startRow = 1, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(wb = wb, file = paste(outPrefix, ".functions.xlsx", sep = ""), overwrite = TRUE)

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
homology_group_heatmap <- function(mat, phy, metadata, hgAn, width, markGenomes){
  
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    # https://github.com/jokergoo/ComplexHeatmap/issues/949
    all(rownames(mat) == phy$tip.label),
    all(rownames(mat) %in% metadata$Genome),
    all(colnames(mat) == hgAn$homology_group_id),
    any(class(phy) == "phylo"),
    length(width) == 2,
    all(is.numeric(width)),
    is.list(markGenomes),
    all(unlist(markGenomes) %in% rownames(mat))
  )
  
  ## homology group heatmap
  ht_hg <- ComplexHeatmap::Heatmap(
    matrix = hgMat,
    name = "hg",
    col = structure(
      viridisLite::viridis(n = min(3, max(hgMat))+1, option = "B"),
      names = seq(0, min(3, max(hgMat)))
    ),
    # heatmap_legend_param = list(legend_height = unit(3, "cm")),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    column_split = hgAn$chr,
    cluster_columns = FALSE, cluster_column_slices = FALSE,
    column_order = hgAn$homology_group_id,
    column_labels = dplyr::pull(hgAn, label, name = homology_group_id),
    show_column_names = TRUE, column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 10),
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE,
    column_title = "Homology groups",
    width = unit(width[2], "cm")
  )
  
  
  # draw(ht_hg)
  
  ## species name key heatmap
  speciesMat <- tibble::tibble(Genome = rownames(mat)) %>% 
    dplyr::left_join(
      y = tibble::enframe(markGenomes) %>% tidyr::unnest(cols = value),
      by = c("Genome" = "value")
    ) %>%
    tidyr::replace_na(replace = list(name = "1")) %>% 
    dplyr::left_join(y = dplyr::select(metadata, Genome, SpeciesName), by = "Genome") %>% 
    tidyr::pivot_wider(
      id_cols = Genome, names_from = SpeciesName,
      values_from = name, values_fill = "0", names_sort = TRUE
    ) %>% 
    tibble::column_to_rownames(var = "Genome") %>% 
    as.matrix()
  
  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == phy$tip.label))
  
  ht_species <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = c("1" = "black", "0" = "white", "compare" = "red", "against" = "green"),
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
# prepare homology group PAV matrix from pan.db
hgMat <- homology_groups_mat(pandb = orgDb, type = "cnv", groups = hgAnnotation$homology_group_id)

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo,
  hgAn = hgAnnotation, width = c(9, 20),
  markGenomes = associatedGenomes
)

pdf(file = paste(outPrefix, ".pheno_hg_association.pdf", sep = ""), width = 14, height = 8)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
)
dev.off()

################################################################################


