#! /usr/bin/env RScript

suppressMessages(library(tidyverse))
suppressMessages(library(ape))

# generate trees for various distance metrics

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
panConf <- confs$data$pangenomes[[pangenome]]

################################################################################
# kmer distance UPGMA tree
df <- suppressMessages(
  readr::read_csv(panConf$db$kmer_classification$KC.100.0$files$mash_dist)
)

distMat <- tibble::column_to_rownames(df, var = "Genomes") %>%
  as.matrix() %>%
  as.dist()

## store UPGMA and NJ trees
# plot(hclust(distMat))
upgmaTree <- ape::as.phylo(hclust(d = distMat, method = "average")) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

ape::write.tree(
  phy = upgmaTree, tree.names = "upgma_tree",
  file = panConf$db$kmer_classification$KC.100.0$files$mash_upgma 
)


njTree <- ape::nj(distMat) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

## set negative length edges => 0
njTree$edge.length[njTree$edge.length < 0] <- 0

ape::write.tree(
  phy = njTree, tree.names = "nj_tree",
  file = panConf$db$kmer_classification$KC.100.0$files$mash_nj
)

################################################################################
