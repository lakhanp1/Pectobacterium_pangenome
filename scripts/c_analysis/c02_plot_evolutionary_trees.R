#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ggstar))
suppressPackageStartupMessages(library(viridisLite))

## A script to visualize phylogenetic trees with some metadata

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser, opt_str = c("-t", "--tree"), type = "character", action = "store",
  help = "tree file in NEWICK format"
)

parser <- optparse::add_option(
  parser, opt_str = c("-n", "--name"), type = "character", action = "store",
  help = "name of the tree"
)

parser <- optparse::add_option(
  parser, opt_str = c("-o", "--out"), type = "character", action = "store",
  help = "output file prefix"
)

parser <- optparse::add_option(
  parser, opt_str = c("-c", "--config"), type = "character", action = "store",
  help = "project config YAML file"
)

opts <- optparse::parse_args(parser)

################################################################################
confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = opts$config)),
  # conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

# ## for test
# confs <- prefix_config_paths(
#   conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
# )
# opts$tree <- confs$data$pangenomes$pectobacterium.v2$db$kmer_classification$KC.100.0

# listviewer::jsonedit(confs)
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
outGroup <- confs$analysis$phylogeny$outgroup

outDir <- confs$analysis$phylogeny$dir

################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## read tree
rawTree <- ape::read.tree(file = opts$tree)

## set negative length edges => 0
rawTree$edge.length[rawTree$edge.length < 0] <- 0

## add data to tree
treeTbl <- as_tibble(rawTree) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata() %>% 
  treeio::root(outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE)

################################################################################
## plotting

pt_tree <- ggtree::ggtree(tr = treeTbl) +
  labs(title = opts$name)

## mark outgroup
pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleName")

## mark species of interest
pt_tree3 <- pt_tree2 +
  ggtree::geom_tiplab(
    mapping = aes(color = SpeciesName, label = nodeLabs),
    size = 3, align = TRUE, linesize = 0.5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(
    values = setNames(
      c("red", "#E57A44", "#088734", "#088734"),
      c(sampleInfoList[[outGroup]]$SpeciesName, "P. brasiliense",
        "P. carotovorum", "P. c. subsp. carotovorum")
    ),
    breaks = NULL,
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()


pt_tree4 <- annotate_ggtree(pt = pt_tree3, offset = 0.25)

# pt_tree4

ggsave(
  plot = pt_tree4, width = 10, height = 20, scale = 2,
  filename = paste(outDir, "/", opts$out, "_tree.pdf", sep = "")
)
