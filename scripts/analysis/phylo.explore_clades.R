#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))

## visualize clades of tree for exploration only

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

phyMethod <- "kmer_nj"
treeConf <- confs$analysis$phylogeny[[phyMethod]]
################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId 
)

## read tree
rawTree <- ape::read.tree(treeConf$files$tree_rooted)

## add data to tree
treeTbl <- as_tibble(rawTree) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "genomeId")) %>%
  treeio::as.treedata()

################################################################################

# pt_tree <- ggtree::ggtree(tr = treeTbl) +
#   labs(title = treeConf$name)
#
# ## mark outgroup
# pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleName")
#
# ## mark species of interest
# pt_tree3 <- pt_tree2 +
#   ggtree::geom_nodelab(
#     mapping = aes(label = node),
#     node = "internal", size = 3, hjust = 1.3
#   ) +
#   ggtree::geom_tiplab(
#     mapping = aes(color = SpeciesName, label = nodeLabs),
#     size = 3, align = TRUE, linesize = 0.5
#   ) +
#   scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
#   scale_color_manual(
#     values = setNames(
#       c("red", "#E57A44", "#088734", "#088734"),
#       c(sampleInfoList[[outGroup]]$SpeciesName, "P. brasiliense",
#         "P. carotovorum", "P. c. subsp. carotovorum")
#     ),
#     breaks = NULL,
#     na.value = "black"
#   ) +
#   ggnewscale::new_scale_color()
#
# pt_tree4 <- annotate_ggtree(pt = pt_tree3, offset = 0.25)
#
#
# ggsave(
#   plot = pt_tree4, width = 10, height = 20, scale = 2,
#   filename = paste(outPrefix, "_tree.pdf", sep = "")
# )

# ggtree::viewClade(tree_view = pt_tree4, node = 626)
# ggtree::scaleClade(tree_view = pt_tree4, node = 626, scale = 8)

subTree2 <- ape::extract.clade(phy = rawTree, node = "n172")
plot.phylo(ladderize(subTree2), show.node.label = T, align.tip.label = T)

ggtree::ggtree(tr = subTree2) +
  ggtree::geom_nodelab(
    mapping = aes(label = label), color = "red", size = 5, hjust = 1.3
  ) +
  ggtree::geom_tiplab(
    mapping = aes(label = label), size = 3, align = TRUE
  )
