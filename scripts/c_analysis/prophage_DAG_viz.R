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

# visualize the prophage network

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
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

sampleInfo <- get_metadata(file = panConf$files$metadata)

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

prophageDf <- suppressMessages(readr::read_tsv(confs$data$genomad$files$prophages)) %>% 
  dplyr::select(prophage_id, prophage_length = length, completeness)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

################################################################################

phageRelations <- suppressMessages(
  readr::read_tsv(file = confs$analysis$prophages$files$prophage_similarity)
) %>% 
  dplyr::mutate(child = stringr::str_split(child, ";")) %>% 
  tidyr::unnest(child) %>% 
  as.data.frame()

parents <- dplyr::filter(phageRelations, !is.na(parent)) %>% 
  dplyr::select(parent) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(nodeType = "root")

nodeType <- dplyr::filter(phageRelations, is.na(parent)) %>% 
  dplyr::select(child) %>% 
  dplyr::left_join(y = parents, by = c("child" = "parent")) %>% 
  tidyr::replace_na(replace = list(nodeType = "single"))

nodes <- dplyr::select(phageRelations, id = child, Genome = childGenome,
                       nHgs = nHgChild, perSharedParent) %>% 
  dplyr::left_join(y = nodeType, by = c("id" = "child")) %>% 
  dplyr::left_join(y = prophageDf, by = c("id" = "prophage_id")) %>% 
  tidyr::replace_na(replace = list(nodeType = "child"))

edges <- dplyr::filter(phageRelations, !is.na(parent)) %>% 
  dplyr::select(
    from = child, to = parent, weight = perSharedParent,
    nSharedHgs, starts_with("perShared"), relation
  )


if(length(setdiff(union(edges$from, edges$to), nodes$id)) != 0) {
  stop(
    "Nodes from edges df are missing in nodes df: ",
    paste(setdiff(union(edges$from, edges$to), nodes$id), collapse = "; ")
  )
}

gf <- igraph::graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)
################################################################################

# there should be only two elements
clique_size_counts(gf)
# largest_cliques(gf)


# components of graph
clu <- igraph::components(gf)
igraph::count_components(gf)
compL <- igraph::groups(clu)

################################################################################
# DAG visualization 
gf %>% 
  plot(vertex.size = 1.5, vertex.label = NA,
       edge.arrow.size = 0.3, edge.color = "black")

gf %>% 
  add_layout_(igraph::as_tree(), component_wise()) %>% 
  plot(vertex.size = 1.5, vertex.label = NA,
       edge.arrow.size = 0.3, edge.color = "black")

gf %>% 
  # add_layout_(igraph::as_star(), igraph::component_wise()) %>% 
  plot(layout = layout_with_dh, vertex.size = 1.5, vertex.label = NA,
       edge.arrow.size = 0.3, edge.color = "black")

# visualization
# layout
lot <- layout_with_mds(graph = gf, dim = 2)
lot <- layout_with_dh(graph = gf)
lot <- layout_with_sugiyama(graph = gf)
lot <- layout_with_lgl(graph = gf)
lot <- layout_with_kk(graph = gf)

plot(gf, layout = lot, vertex.size = 1.5, vertex.label = NA,
     edge.arrow.size = 0.3, edge.color = "black")

# each component layout first and then arranged
lot <- igraph::layout_components(graph = gf, layout = layout_with_lgl)
plot(gf, layout = lot, vertex.size = 1.5, vertex.label = NA,
     edge.arrow.size = 0.3, edge.color = "black")


(pt_net <- ggraph(gf, layout = 'kk') +
    geom_edge_link(
      # mapping = aes(color = weight),
      arrow = arrow(length = unit(1, 'mm')),
      end_cap = circle(1, 'mm')
    ) +
    # scale_colour_viridis_c(option = "A") +
    # ggnewscale::new_scale_colour() +
    geom_node_point(mapping = aes(color = nodeType, size = nodeType), alpha = 0.7) +
    scale_color_manual(
      values = c("child" = "black", "root" = "#d95f02", single = "#7570b3")
    ) +
    scale_size_manual(
      values = c("child" = 1.5, "root" = 3, single = 1.5)
    ) +
    labs(title = "pangenome prophage hierarchy") +
    theme_void(base_size = 18) +
    theme(
      legend.position = "right"
    ))

ggsave(filename = paste(outDir, "/prophage_clusters.pdf", sep = ""),
       width = 12, height = 10)

################################################################################

# prophage homology groups on pangenome phylogeny
proReps <- dplyr::filter(nodes, nodeType != "child") %>% 
  dplyr::select(prophage_id = id, nodeType) %>% 
  dplyr::left_join(proHgs, by = "prophage_id") %>% 
  dplyr::select(-contig_id)

# function to get homology group counts in each genome of a pangenome
homology_groups_genome_counts <- function(h, pandb){
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
)


countMat <- tibble::tibble(genome = rawTree$tip.label) %>% 
  dplyr::left_join(y = proHgCounts, by = c("genome")) %>% 
  tibble::column_to_rownames(var = "genome") %>% 
  as.matrix()

hist(matrixStats::colMaxs(countMat))
hist(countMat)

speciesMat <- tibble::tibble(Genome = rownames(countMat)) %>% 
  dplyr::left_join(y = dplyr::select(sampleInfo, Genome, SpeciesName), by = "Genome") %>% 
  dplyr::mutate(species = 1) %>% 
  tidyr::pivot_wider(
    id_cols = Genome, names_from = SpeciesName,
    values_from = species, values_fill = 0, names_sort = TRUE
  ) %>% 
  tibble::column_to_rownames(var = "Genome") %>% 
  as.matrix()

prophageAnDf <- tibble::tibble(prophage_id = colnames(countMat)) %>% 
  dplyr::left_join(y = nodes, by = c("prophage_id" = "id"))

## ensure the row order is same: this is because of a bug in ComplexHeatmap
stopifnot(
  all(rownames(countMat) == rawTree$tip.label),
  all(rownames(speciesMat) == rawTree$tip.label)
)

# top annotation
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
  cluster_columns = TRUE,
  show_column_names = FALSE,
  show_row_dend = TRUE, show_column_dend = FALSE,
  row_dend_width = unit(3, "cm"),
  show_row_names = FALSE,
  column_title = "Prophage homology groups / genome",
  width = unit(20, "cm")
)


# species key heatmap
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
  width = unit(8, "cm")
)

htList <- ht_species + ht_count

pdf(file = paste(outDir, "/prophage_hg_summary.pdf", sep = ""), width = 15, height = 9)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg_count",
  row_dend_side = "left",
  merge_legends = TRUE
)
dev.off()



