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

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser,
  opt_str = c("-t", "--tree"), type = "character", action = "store",
  help = "tree file in NEWICK format"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-n", "--name"), type = "character", action = "store",
  help = "name of the tree from YAML phylogeny tags"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-c", "--config"), type = "character", action = "store",
  help = "project config YAML file"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("--save_leaf_order"), type = "logical", action = "store_true",
  help = "LOGICAL: Save the order of tip lables from phylogenetic tree"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$name), is.null(opts$config), is.null(opts$tree))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# ## for test
# opts$config <- "project_config.yaml"
# opts$tree <- "./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/kmer_classification.100.0/kmer.upgma.newick"
# opts$name <- "kmer_upgma"
# opts$save_leaf_order <- TRUE
# opts$tree <- "./data/pangenomes/pectobacterium.v2/pectobacterium.v2.DB/core_snp_tree/informative.fasta.treefile"
# opts$name <- "core_snp_ml"
# #######

################################################################################

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = opts$config)),
  dir = "."
)

stopifnot(
  file.exists(opts$tree),
  !is.null(opts$name),
  has_name(confs$analysis$phylogeny, opts$name),
  !is.null(confs$analysis$phylogeny[[opts$name]]$files$tree)
)

# listviewer::jsonedit(confs)
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
outGroup <- confs$analysis$phylogeny$outgroup

outDir <- confs$analysis$phylogeny[[opts$name]]$dir

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## read tree
rawTree <- ape::read.tree(file = opts$tree)

## set negative length edges => 0
rawTree$edge.length[rawTree$edge.length < 0] <- 0

rawTree <- ape::ladderize(rawTree) %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

ape::write.tree(
  phy = rawTree, tree.names = opts$name,
  file = confs$analysis$phylogeny[[opts$name]]$files$tree
)


rootedTr <- ape::root(
  phy = rawTree, outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE
) %>%
  ape::ladderize()

ape::write.tree(
  phy = rootedTr, tree.names = opts$name,
  file = confs$analysis$phylogeny[[opts$name]]$files$tree_rooted
)

## add data to tree
treeTbl <- as_tibble(rootedTr) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata()

################################################################################
## plotting

pt_tree <- ggtree::ggtree(tr = treeTbl) +
  labs(title = opts$name)

## mark outgroup
pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleName")

## mark species of interest
pt_tree3 <- pt_tree2 +
  ggtree::geom_nodelab(
    mapping = aes(label = label),
    node = "internal", size = 3, hjust = 1.3, color = "red"
  ) +
  ggtree::geom_tiplab(
    mapping = aes(color = SpeciesName, label = nodeLabs),
    size = 3, align = TRUE, linesize = 0.5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(
    values = setNames(
      c("red", "#E57A44", "#088734", "#088734"),
      c(
        sampleInfoList[[outGroup]]$SpeciesName, "P. brasiliense",
        "P. carotovorum", "P. c. subsp. carotovorum"
      )
    ),
    breaks = NULL,
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()


pt_tree4 <- annotate_ggtree(pt = pt_tree3, offset = 0.25)

# pt_tree4

ggsave(
  plot = pt_tree4, width = 10, height = 20, scale = 2,
  filename = paste(outDir, "/", opts$name, ".rooted_tree.pdf", sep = "")
)

################################################################################
# optionally, save species order for plotting species key heatmap
if(opts$save_leaf_order){
  leafOrder <- dplyr::arrange(.data = pt_tree$data, y) %>%
    dplyr::filter(isTip)
  
  speciesOrder <- dplyr::select(leafOrder, SpeciesName, y, type_material) %>%
    dplyr::group_by(SpeciesName) %>%
    dplyr::arrange(type_material, y, .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(y)) %>%
    # dplyr::mutate(SpeciesName = forcats::as_factor(SpeciesName)) %>%
    dplyr::pull(SpeciesName)
  
  readr::write_tsv(
    x = tibble::tibble(SpeciesName = speciesOrder),
    file = confs$analysis$phylogeny[[opts$name]]$files$species_order
  )
  
  # add nodepath to the pangenome metadata file
  nodePathCol <- paste('nodepath.', opts$name, sep = "")
  nodePaths <- nodepath_df(phy = rootedTr) %>% 
    dplyr::rename(
      Genome = tip,
      !!nodePathCol := nodepath
    ) %>% 
    dplyr::mutate(Genome = as.numeric(Genome))
  
  panMeta <- suppressMessages(readr::read_csv(panConf$files$metadata))
  
  if(tibble::has_name(panMeta, nodePathCol)){
    panMeta %<>% dplyr::select(-!!nodePathCol) 
  }
  
  panMeta %<>% dplyr::left_join(nodePaths, by = "Genome")
  
  readr::write_csv(
    x = panMeta, file = panConf$files$metadata
  )
  
  ## write metadata to excel
  wb <- openxlsx::createWorkbook()
  currentSheet <- "metadata"
  openxlsx::addWorksheet(wb, sheetName = currentSheet)
  openxlsx::writeDataTable(
    wb = wb, sheet = currentSheet, withFilter = TRUE, keepNA = TRUE,
    na.string = "NA",
    x = panMeta
  )
  openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)
  
  # openxlsx::openXL(wb)
  openxlsx::saveWorkbook(
    wb = wb,
    file = file.path(panConf$dir, "pangenome_metadata.xlsx"), overwrite = TRUE
  )
}

################################################################################




