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

# parser <- optparse::add_option(
#   parser, default = FALSE,
#   opt_str = c("-r", "--root"), type = "logical", action = "store_true",
#   help = "LOGICAL: whether to root the tree"
# )

parser <- optparse::add_option(
  parser, default = FALSE,
  opt_str = c("-b", "--bootstrap"), type = "logical", action = "store_true",
  help = "LOGICAL: whether to parse the bootstrap values"
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
  parser, default = FALSE,
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
# opts$bootstrap <- TRUE
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

sampleInfo <- get_metadata(
  file = confs$data$pangenomes[[pangenome]]$files$metadata,
  genus = confs$genus
) %>%
  dplyr::mutate(
    genomeString = paste(SpeciesName, " | (", strain, ")", sep = "")
  )


sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

# read tree
if(opts$bootstrap){
  rawTree <- treeio::read.newick(file = opts$tree, node.label = "support")
  bootstrapData <- dplyr::rename(rawTree@data, bootstrap = support)
  treePhy <- rawTree@phylo

} else{
  treePhy <- ape::read.tree(file = opts$tree)
}

# add a 'g_' prefix if tree tips are numeric
if(any(grepl(pattern = "^\\d+$", treePhy$tip.label))){
  treePhy$tip.label <- paste("g_", treePhy$tip.label, sep = "")
}

## set negative length edges => 0
treePhy$edge.length[treePhy$edge.length < 0] <- 0

treePhy <- ape::ladderize(treePhy) %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

rootedTr <- ape::root(
  phy = treePhy, outgroup = sampleInfoList[[outGroup]]$genomeId,
  edgelabel = TRUE, resolve.root = TRUE
) %>%
  ape::ladderize() %>%
  ape::makeNodeLabel(method = "number", prefix = "n")

# save both trees in BEAST format
treePhyData <- tidytree::as.treedata(treePhy)
rootedTreeData <- tidytree::as.treedata(rootedTr)

if(opts$bootstrap){
  treePhyData <- treeio::full_join(
    x = treePhy, y = bootstrapData, by = "node"
  )

  rootedTreeData <- treeio::full_join(
    x = rootedTr, y = bootstrapData, by = "node"
  )

}

# store tree in BEAST format
treeio::write.beast(
  treedata = treePhyData, tree.name = opts$name,
  file = confs$analysis$phylogeny[[opts$name]]$files$tree
)

treeio::write.beast(
  treedata = rootedTreeData, tree.name = opts$name,
  file = confs$analysis$phylogeny[[opts$name]]$files$tree_rooted
)

## add data to tree
treeTbl <- dplyr::full_join(
  x = rootedTr, y = sampleInfo, by = c("label" = "genomeId")
)

################################################################################
## plotting

pt_tree <- ggtree::ggtree(tr = treeTbl) +
  geom_treescale(
    x = 0, y = nrow(sampleInfo) * 0.8,
    fontsize = 8, linesize = 2, offset = 4
  ) +
  labs(title = opts$name)

## mark outgroup
pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleId")

## mark species of interest
pt_tree3 <- pt_tree2 +
  ggtree::geom_nodelab(
    mapping = aes(label = label),
    node = "internal", size = 3, hjust = 1.3, color = "red"
  ) +
  ggtree::geom_tiplab(
    mapping = aes(label = label, color = type_material),
    size = 3, align = TRUE, linesize = 0.5
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(
    values = c("type strain" = "blue"),
    name = "Type strain",
    guide = guide_legend(order = 2),
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()

## add metadata
pt_tree4 <- pt_tree3 +
  ## SpeciesName
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, label = genomeString),
    geom = "geom_text", size = 3, hjust = "left",
    pwidth = 0.6, offset = 0.06
  ) +
  ggnewscale::new_scale_color() +
  ## virulence PCR result
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "vir_pcr", color = virulence_pcr),
    geom = "geom_point",
    pwidth = 0.01, offset = 0.2
  ) +
  scale_color_manual(
    name = "Blackleg PCR",
    values = c("positive" = "red", "negative" = "green"),
    labels = c("positive" = "BL +ve", "negative" = "BL -ve"),
    guide = guide_legend(order = 3),
    na.value = alpha("white", 0)
  ) +
  ggnewscale::new_scale_color() +
  ## virulence phenotype
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "virulence", color = virulence),
    geom = "geom_point", shape = 17, size = 2,
    pwidth = 0.01, offset = 0.02
  ) +
  scale_color_manual(
    name = "Blackleg phenotype",
    values = c("virulent" = "red", "avirulent" = "green"),
    labels = c("virulent" = "BL-causing", "avirulent" = "BL-non-causing"),
    guide = guide_legend(order = 4),
    na.value = alpha("white", 0)
  ) +
  ## collection year
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, label = collection_year),
    geom = "geom_text", size = 3, hjust = "left",
    pwidth = 0.2
  ) +
  ## country
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, label = geo_loc_country),
    geom = "geom_text", size = 3, hjust = "left",
    pwidth = 0.2, offset = 0.05
  )  +
  ggnewscale::new_scale_color() +
  ## BioSampleAccn
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, label = BioSampleAccn, color = source),
    geom = "geom_text", size = 3, hjust = "left",
    pwidth = 0.2, offset = 0.15
  ) +
  scale_color_manual(
    values = c("NCBI" = "black", "NAK" = "#377eb8", "NIVIP" = "#ff7f00"),
    labels = c("NCBI" = "NCBI", "NAK" = "NAK (this study)", "NIVIP" = "NIVIP (this study)"),
    guide = guide_legend(order = 5)
  )


# pt_tree4

ggsave(
  plot = pt_tree4, width = 10, height = 20, scale = 2,
  filename = paste(outDir, "/", opts$name, ".rooted_tree.pdf", sep = "")
)

sampleInfo %<>% dplyr::select(-genomeString)

################################################################################
# optionally, save species order for plotting species key heatmap
if(opts$save_leaf_order){

  leafOrder <- tibble:::enframe(
    ggtree::get_taxa_name(pt_tree), name = "tipOrder", value = "genomeId"
  ) %>%
    dplyr::left_join(y = sampleInfo, by = "genomeId")


  speciesOrder <- dplyr::select(leafOrder, SpeciesName, tipOrder, type_material) %>%
    dplyr::group_by(SpeciesName) %>%
    dplyr::arrange(type_material, tipOrder, .by_group = TRUE) %>%
    dplyr::slice(1L) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(desc(tipOrder)) %>%
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
      genomeId = tip,
      !!nodePathCol := nodepath
    )

  if(tibble::has_name(sampleInfo, nodePathCol)){
    sampleInfo %<>% dplyr::select(-!!nodePathCol)
  }

  sampleInfo %<>% dplyr::left_join(nodePaths, by = "genomeId") %>%
    dplyr::select(Genome, genomeId, everything(), -nodeLabs)

  readr::write_csv(
    x = sampleInfo, file = panConf$files$metadata
  )

  ## write metadata to excel
  wb <- openxlsx::createWorkbook()
  currentSheet <- "metadata"
  openxlsx::addWorksheet(wb, sheetName = currentSheet)
  openxlsx::writeDataTable(
    wb = wb, sheet = currentSheet, withFilter = TRUE, keepNA = TRUE,
    na.string = "NA",
    x = sampleInfo
  )
  openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

  # openxlsx::openXL(wb)
  openxlsx::saveWorkbook(
    wb = wb,
    file = file.path(panConf$dir, "pangenome_metadata.xlsx"), overwrite = TRUE
  )
}

################################################################################




