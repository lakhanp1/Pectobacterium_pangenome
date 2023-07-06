#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))

# post-processing of prophage DAG to add additional information:
# subgraph, root node for each subgraph

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

treeMethod <- "kmer_nj"     #ani_upgma, kmer_nj

outDir <- confs$analysis$prophages$dir

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

# get root-leaf node path as string for filtering purpose
rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)
nodePaths <- nodepath_df(phy = rawTree) %>% 
  dplyr::rename(Genome = tip) %>% 
  dplyr::mutate(Genome = as.numeric(Genome))

prophageDf <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>%
  dplyr::select(prophage_id,
    prophage_length = length, prophage_taxonomy = taxonomy,
    completeness, checkv_quality, sampleId, SpeciesName
  )

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>% 
  dplyr::select(prophage_id, hgs)

# proHgL <- purrr::transpose(proHgs) %>%
#   purrr::set_names(nm = purrr::map(., "prophage_id"))

sharedStats <- suppressMessages(readr::read_tsv(
  file = confs$analysis$prophages$files$hg_similarity
))

################################################################################
# store network format table
phageRelations <- dplyr::mutate(sharedStats, child = stringr::str_split(child, ";")) %>%
  tidyr::unnest(child)

parents <- dplyr::filter(phageRelations, !is.na(parent)) %>%
  dplyr::select(parent) %>%
  dplyr::distinct() %>%
  dplyr::mutate(nodeType = "root")

nodeType <- dplyr::filter(phageRelations, is.na(parent)) %>%
  dplyr::select(child) %>%
  dplyr::left_join(y = parents, by = c("child" = "parent")) %>%
  tidyr::replace_na(replace = list(nodeType = "singleton"))

phageRelations %<>%
  dplyr::left_join(y = nodeType, by = c("child")) %>%
  tidyr::replace_na(replace = list(nodeType = "child")) %>%
  dplyr::mutate(
    nodeType = forcats::fct_relevel(.f = nodeType, "root", "singleton")
  ) %>%
  dplyr::rename(prophage_id = child, Genome = childGenome, nHgs = nHgChild) %>%
  dplyr::relocate(nodeType, .after = prophage_id) %>% 
  dplyr::left_join(y = nodePaths, by = "Genome") %>% 
  dplyr::left_join(y = proHgs, by = "prophage_id")


# make igraph
nodes <- dplyr::select(
  phageRelations, prophage_id, nodeType, Genome,
  nHgs, perSharedParent
) %>%
  dplyr::left_join(y = prophageDf, by = "prophage_id") %>%
  dplyr::rename(id = prophage_id)

edges <- dplyr::filter(phageRelations, !is.na(parent)) %>%
  dplyr::select(
    from = prophage_id, to = parent, perSharedParent,
    nSharedHgs, perSharedChild, relation
  ) %>%
  dplyr::mutate(weight = perSharedParent)


if (length(setdiff(union(edges$from, edges$to), nodes$id)) != 0) {
  stop(
    "Nodes from edges df are missing in nodes df: ",
    paste(setdiff(union(edges$from, edges$to), nodes$id), collapse = "; ")
  )
}

gf <- igraph::graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)

# components of graph
igraph::count_components(gf)
clu <- igraph::components(gf)
compDf <- tibble::enframe(clu$membership, name = "prophage_id", value = "graph_component") %>%
  dplyr::group_by() %>%
  dplyr::add_count(graph_component, name = "component_size")

# identify root prophage for each component
phageRelations <- dplyr::left_join(
  phageRelations, compDf,
  by = "prophage_id"
)

compRoots <- dplyr::filter(phageRelations, nodeType != "child") %>%
  dplyr::select(
    root_id = prophage_id, root_nHgs = nHgs, root_hgs = hgs, graph_component
  ) 

networkDf <- dplyr::left_join(
  x = phageRelations, y = compRoots, by = "graph_component"
) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    root_shared = length(intersect(hgs, root_hgs)),
    root_jaccard = root_shared / length(union(hgs, root_hgs)),
    root_contentDissimilarity = 1 - (mean(c(root_shared/nHgs, root_shared/root_nHgs))),
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(-c(hgs, root_hgs)) %>% 
  dplyr::arrange(graph_component, nodeType)

readr::write_tsv(
  x = networkDf,
  file = confs$analysis$prophages$files$network
)
