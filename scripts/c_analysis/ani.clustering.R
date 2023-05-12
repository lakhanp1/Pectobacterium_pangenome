suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

## generate UPGMA and NJ tree using ANI for pangenome

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
outGroup <- confs$analysis$phylogeny$outgroup

outDir <- confs$analysis$phylogeny$ani$dir

if(!dir.exists(outDir)) dir.create(outDir)

pt_theme <- theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, face = "bold")
  )

################################################################################
sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

genomeIds <- dplyr::pull(sampleInfo, Genome, name = sampleId)

################################################################################
## process ANI data for pangenome and store ANI and ANI-distance matrices
aniDf <- suppressMessages(readr::read_tsv(
  file = confs$analysis$ANI$files$fastani_out,
  col_names = c("id1", "id2", "ani", "mapped", "total")
))

aniDf %<>% dplyr::mutate(
  dplyr::across(
    .cols = c(id1, id2),
    .fns = ~stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1"),
  ),
  dist = 1 - (ani/100)
) %>% 
  dplyr::mutate(
    g1 = genomeIds[id1],
    g2 = genomeIds[id2]
  ) %>% 
  dplyr::filter(!is.na(g1) & !is.na(g2)) %>% 
  dplyr::arrange(g1, g2)

aniDist <- tidyr::pivot_wider(
  data = aniDf,
  id_cols = "g1",
  names_from = "g2",
  values_from = "dist"
)

readr::write_tsv(
  x = aniDist, file = confs$analysis$phylogeny$ani$files$ani_distance
)

distMat <- tibble::column_to_rownames(aniDist, var = "g1") %>% 
  as.matrix() %>% 
  as.dist()

aniMat <- tidyr::pivot_wider(
  data = aniDf,
  id_cols = "g1",
  names_from = "g2",
  values_from = "ani"
)

readr::write_tsv(
  x = aniMat, file = confs$analysis$phylogeny$ani$files$ani_matrix
)

aniMat <- tibble::column_to_rownames(aniMat, var = "g1") %>% 
  as.matrix()

if(!all(rownames(as.matrix(distMat)) == rownames(aniMat))){
  stop("rownames did not match")
}
################################################################################
## store UPGMA and NJ trees

# plot(hclust(distMat))
treeUpgma <- ape::as.phylo(hclust(d = distMat, method = "average")) %>% 
  ape::ladderize() %>% 
  ape::makeNodeLabel(method = "number", prefix = "n")

ape::write.tree(
  phy = treeUpgma, tree.names = "ani_upgma",
  file = confs$analysis$phylogeny$ani_upgma$files$tree
)

# plot(ape::root(phy = treeUpgma, outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE))
# nodelabels()

treeNj <- ape::nj(distMat) %>%
  ape::ladderize() %>% 
  ape::makeNodeLabel(method = "number", prefix = "n")

## set negative length edges => 0
treeNj$edge.length[treeNj$edge.length < 0] <- 0

ape::write.tree(
  phy = treeNj, tree.names = "ANI_NJ",
  file = confs$analysis$phylogeny$ani_nj$files$tree
)



################################################################################
## add data to tree
treeTbl <- as_tibble(treeUpgma) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata() %>% 
  treeio::root(outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE)

pt_treeUpgma <- ggtree::ggtree(
  tr = treeTbl
)

## get species order to arrange the species key columns
leafOrder <- dplyr::arrange(.data = pt_treeUpgma$data, y) %>% 
  dplyr::filter(isTip)

speciesOrder <- dplyr::select(leafOrder, SpeciesName, y, type_material) %>% 
  dplyr::group_by(SpeciesName) %>% 
  dplyr::arrange(type_material, y, .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% dplyr::arrange(desc(y)) %>% 
  # dplyr::mutate(SpeciesName = forcats::as_factor(SpeciesName)) %>% 
  dplyr::pull(SpeciesName)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder)
)

readr::write_tsv(
  x = tibble::tibble(SpeciesName = speciesOrder),
  file = confs$analysis$phylogeny$ani_upgma$files$species_order
)

################################################################################
## visualize ANI data: all tree
#' Plot ANI heatmap with provided clustering
#'
#' @param mat ANI score matrix
#' @param phy A phylo object
#' @param metadata metadata for genomes. Must have SpeciesName and Genome column
#' @param width Vector of length 2 for heatmap widths
#'
#' @return ComplexHeatmap with layout: species key heatmap | ANI heatmap
#' @export
#'
#' @examples
plot_species_ANI_heatmap <- function(mat, phy, metadata, width){
  
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    all(rownames(mat) == phy$tip.label),
    setequal(rownames(mat), colnames(mat)),
    all(rownames(mat) %in% metadata$Genome),
    any(class(phy) == "phylo"),
    length(width) == 2,
    all(is.numeric(width))
  )
  
  ## ANI heatmap
  ht_ani <- ComplexHeatmap::Heatmap(
    matrix = mat,
    name = "ani",
    col = circlize::colorRamp2(
      breaks = c(80, 85, 90, 91, 92, 92.5, 93, 93.5, 94, 95, 96, 97, 99),
      # breaks = c(80, 85, 90, 91, 91.5, 92, 92.5, 93, 93.5, 94, 96, 97, 99),
      colors = viridisLite::viridis(n = 13, option = "B")
    ),
    heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm")),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    cluster_columns = ape::as.hclust.phylo(phy), column_dend_reorder = FALSE,
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE, show_column_names = FALSE,
    column_title = "ANI scores",
    width = unit(width[2], "cm")
  )
  
  # draw(ht_ani)
  ## species name key heatmap
  speciesMat <- get_species_key_data(genomes = rownames(mat), metadata = metadata, type = "wide")
  
  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == phy$tip.label))
  
  ht_species <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = c("1" = "black", "0" = "white"),
    cluster_columns = FALSE,
    column_split = 1:ncol(speciesMat), cluster_column_slices = FALSE,
    border = TRUE, column_gap = unit(0, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_title = "Species key",
    width = unit(width[1], "cm")
  )
  
  # draw(ht_species)
  
  htList <- ht_species + ht_ani
  
  return(htList)
  
}

################################################################################
## ANI heatmap for all genomes
htList <- plot_species_ANI_heatmap(
  mat = aniMat, phy = treeUpgma, metadata = sampleInfo, width = c(12,18)
)

png(filename = file.path(outDir, "ANI_heatmap.png"), width = 5000, height = 2800, res = 350)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "ani",
  row_dend_side = "left"
)
dev.off()

################################################################################
## ANI heatmap for subset of genomes

## ANI heatmap for inhouse strains
# inhouseNodes <- dplyr::filter(sampleInfo, source %in% c("NAK", "NVWA"))
# subsetName <- "ANI_inhouse"

## ANI heatmap for P. brasiliense strains
inhouseNodes <- dplyr::filter(sampleInfo, SpeciesName == "P. brasiliense")
subsetName <- "ANI_PBrasiliense"

subTree2 <- ape::keep.tip(phy = treeUpgma, tip = inhouseNodes$Genome)
subAni2 <- aniMat[subTree2$tip.label, subTree2$tip.label] %>% 
  tibble::as_tibble(rownames = "g1") %>% 
  tidyr::pivot_longer(
    cols = -g1,
    names_to = "g2", values_to = "ANI"
  )

inhouseTreeTbl <- as_tibble(subTree2) %>%
  dplyr::full_join(y = inhouseNodes, by = c("label" = "Genome")) %>%
  treeio::as.treedata()

## generate ggtree plot
pt_inhouseTree <- ggtree::ggtree(
  tr = inhouseTreeTbl
) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ## virulence phenotype
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "virulence", color = virulence),
    geom = "geom_point", shape = 17, size = 2,
    pwidth = 0.01, offset = 0.1
  ) +
  scale_color_manual(
    values = c("virulent" = "red", "avirulent" = "green"),
    na.value = alpha("white", 0)
  ) +
  ggnewscale::new_scale_color() +
  ## virulence PCR result
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "vir_pcr", color = virulence_pcr),
    geom = "geom_point",
    pwidth = 0.01, offset = 0.1
  ) +
  scale_color_manual(
    values = c("positive" = "red", "negative" = "green"),
    na.value = alpha("white", 0)
  )

# # mark specific nodes: SpeciesName
# pt_inhouseTree <- pt_inhouseTree +
#   ggnewscale::new_scale_color() +
#   ggtree::geom_tippoint(
#     mapping = aes(subset = c(SpeciesName == "P. brasiliense")),
#     color = "blue"
#   ) 

# mark specific nodes: source
pt_inhouseTree <- pt_inhouseTree +
  ggnewscale::new_scale_color() +
  ggtree::geom_tippoint(
    mapping = aes(subset = c(source %in% c("NAK", "NVWA"))),
    color = "blue"
  )


# species name key heatmap
speciesKyeDf <- get_species_key_data(
  genomes = inhouseNodes$Genome, metadata = sampleInfo, type = "long"
)

pt_spKey <- species_key_plot(keyDf = speciesKyeDf)

# ANI heatmap
subAni2 <- dplyr::mutate(
  subAni2,
  g1 = forcats::fct_relevel(g1, ggtree::get_taxa_name(pt_inhouseTree)),
  g2 = forcats::fct_relevel(g2, ggtree::get_taxa_name(pt_inhouseTree))
)
pt_ani <- ani_heatmap(aniDf = subAni2)

## arrange plots one by one
pt_all <- pt_spKey %>% aplot::insert_left(pt_inhouseTree, width = 0.5) %>%
  aplot::insert_right(pt_ani, width = 2)


ggsave(
  plot = pt_all, width = 18, height = 10,
  filename = paste(outDir, "/", subsetName, ".heatmap.pdf", sep = "")
)


################################################################################

