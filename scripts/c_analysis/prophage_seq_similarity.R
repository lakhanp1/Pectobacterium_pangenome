#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(corrplot))

# summarize prophage ANI and MASH results
# cluster prophages based on ANI distance

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/genome_comparison.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- confs$analysis$prophages$dir

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

prophageDf <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>%
  dplyr::select(-Genome)

phageRelations <- suppressMessages(
  readr::read_tsv(file = confs$analysis$prophages$files$network)
) %>%
  dplyr::select(id = prophage_id, nodeType, nHgs) %>%
  dplyr::left_join(y = prophageDf, by = c("id" = "prophage_id")) %>%
  dplyr::mutate(
    nodeLabs = paste(
      id, "|", "| hgs =", nHgs, "| length =", length,
      "|", SpeciesName
    )
  ) %>%
  dplyr::filter(nodeType != "child")

nodeOfInterest <- dplyr::filter(phageRelations, nodeType != "child") %>%
  dplyr::select(id)

################################################################################
# process ANI data for phages and store ANI and ANI-distance matrices
aniDf <- suppressMessages(readr::read_tsv(
  file = confs$data$prophages$files$ani,
  col_names = c("id1", "id2", "ani", "mapped", "total")
)) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(id1, id2),
      .fns = ~ stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1"),
    ),
    dist = 1 - (ani / 100)
  )

aniDist <- dplyr::left_join(
  x = nodeOfInterest,
  y = tidyr::pivot_wider(
    data = aniDf,
    id_cols = "id1",
    names_from = "id2",
    values_from = "dist",
    values_fill = 1
  ),
  by = c("id" = "id1")
) %>%
  dplyr::select(id, all_of(nodeOfInterest$id))

readr::write_tsv(
  x = aniDist, file = confs$analysis$prophages$files$ani_distance
)

aniDistMat <- tibble::column_to_rownames(aniDist, var = "id") %>%
  as.matrix() %>%
  as.dist()

## store UPGMA and NJ trees
# plot(hclust(distMat))
aniUpgma <- ape::as.phylo(hclust(d = aniDistMat, method = "average")) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

ape::write.tree(
  phy = aniUpgma, tree.names = "ani_upgma",
  file = confs$analysis$prophages$files$ani_upgma
)

# plot(ape::root(phy = treeUpgma, outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE))
# nodelabels()

aniNj <- ape::nj(aniDistMat) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

## set negative length edges => 0
aniNj$edge.length[aniNj$edge.length < 0] <- 0

ape::write.tree(
  phy = aniNj, tree.names = "ANI_NJ",
  file = confs$analysis$prophages$files$ani_nj
)

################################################################################
# process MASH data and generate trees
mashDf <- suppressMessages(readr::read_tsv(
  file = confs$data$prophages$files$mash,
  col_names = c("id1", "id2", "dist", "pvalue", "mapped")
)) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(id1, id2),
      .fns = ~ stringr::str_replace(
        string = .x, pattern = ".*/(.*).fna", replacement = "\\1"
      ),
    )
  )

# mashDf %>%
#   dplyr::filter(id1 %in% !!nodeOfInterest$id, id2 %in% !!nodeOfInterest$id) %>%
#   dplyr::filter(id1 != id2) %>%
#   dplyr::arrange(dist, id1, id2) %>%
#   dplyr::left_join(
#     y = dplyr::select(prophageDf, id1 = prophage_id, species1 = SpeciesName),
#     by = "id1"
#   ) %>%
#   dplyr::left_join(
#     y = dplyr::select(prophageDf, id2 = prophage_id, species2 = SpeciesName),
#     by = "id2"
#   ) %>%
#   clipr::write_clip()

mashDist <- dplyr::left_join(
  x = nodeOfInterest,
  y = tidyr::pivot_wider(
    data = mashDf,
    id_cols = "id1",
    names_from = "id2",
    values_from = "dist",
    values_fill = 1
  ),
  by = c("id" = "id1")
) %>%
  dplyr::select(id, all_of(nodeOfInterest$id))


readr::write_tsv(
  x = mashDist, file = confs$analysis$prophages$files$mash_distance
)

mashDistMat <- tibble::column_to_rownames(mashDist, var = "id") %>%
  as.matrix() %>%
  as.dist()

## store UPGMA and NJ trees
# plot(hclust(distMat))
mashUpgma <- ape::as.phylo(hclust(d = mashDistMat, method = "average")) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

ape::write.tree(
  phy = mashUpgma, tree.names = "mash_upgma",
  file = confs$analysis$prophages$files$mash_upgma
)


mashNj <- ape::nj(mashDistMat) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

## set negative length edges => 0
mashNj$edge.length[mashNj$edge.length < 0] <- 0

ape::write.tree(
  phy = mashNj, tree.names = "mash_nj",
  file = confs$analysis$prophages$files$mash_nj
)

################################################################################
# correlate trees
## cophenetic distance calculation for NJ and UPGMA tree
copheneticDist <- tibble::tibble(
  aniDist = as.vector(aniDistMat),
  aniUpgma = as.vector(as.dist(cophenetic(aniUpgma))),
  aniNj = as.vector(as.dist(cophenetic(aniNj))),
  mashDist = as.vector(mashDistMat),
  mashUpgma = as.vector(as.dist(cophenetic(mashUpgma))),
  mashNj = as.vector(as.dist(cophenetic(mashNj)))
)

M <- cor(copheneticDist)

corrplot::corrplot(
  M,
  type = "lower", tl.col = "black",
  cl.ratio = 0.2, tl.srt = 45, col = COL2("BrBG"), addCoef.col = "white"
)


################################################################################
## add data to tree
treeTbl <- treeio::as_tibble(mashUpgma) %>%
  dplyr::left_join(y = phageRelations, by = c("label" = "id")) %>%
  treeio::as.treedata()

pt_treeUpgma <- ggtree::ggtree(
  tr = treeTbl
) +
  ggtree::geom_tiplab(
    mapping = aes(color = SpeciesName, label = nodeLabs),
    size = 2, align = TRUE, linesize = 0.5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.5))) +
  scale_color_manual(
    values = setNames(
      viridisLite::viridis(n = 5),
      c(
        "P. brasiliense", "P. carotovorum", "P. c. subsp. carotovorum",
        "P. versatile", "P. parmentieri"
      )
    ),
    breaks = NULL,
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()

ggsave(
  plot = pt_treeUpgma, width = 8, height = 20, scale = 2,
  filename = paste(outDir, "/prophages.MASH.UPGMA_tree.pdf", sep = "")
)
################################################################################

# ANI
ht_ani <- plot_species_ANI_heatmap(
  mat = 100 * (1 - as.matrix(aniDistMat)), phy = aniUpgma, speciesInfo = NULL,
  col = circlize::colorRamp2(
    breaks = c(0, 50, 80, 90, 91, 92, 93, 93.5, 94, 95, 96, 97, 99),
    colors = viridisLite::viridis(n = 13, option = "B")
  ),
  name = "ani"
)

pdf(
  file = file.path(outDir, "prophage_ANI_heatmap.pdf"),
  width = 10, height = 10
)
ComplexHeatmap::draw(
  object = ht_ani,
  main_heatmap = "ani",
  row_dend_side = "left",
  heatmap_legend_side = "bottom"
)
dev.off()


# MASH
ht_mash <- plot_species_ANI_heatmap(
  mat = as.matrix(mashDistMat), phy = mashUpgma, speciesInfo = NULL,
  name = "mash",
  col = circlize::colorRamp2(
    breaks = c(0, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1),
    # breaks = seq(0, 1, length.out = 13),
    colors = viridisLite::viridis(n = 13, option = "B")
  ),
)

pdf(
  file = file.path(outDir, "prophage_MASH_heatmap.pdf"),
  width = 10, height = 10
)
ComplexHeatmap::draw(
  object = ht_mash,
  main_heatmap = "mash",
  row_dend_side = "left",
  heatmap_legend_side = "bottom"
)
dev.off()

################################################################################
