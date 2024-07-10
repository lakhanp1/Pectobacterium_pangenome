#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

## visualize homology groups PAV for given homology group sets

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/heatmap_utils.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser,
  opt_str = c("-t", "--tree"), type = "character", action = "store",
  help = "species tree name. Tree file will be fetched from the config file."
)

parser <- optparse::add_option(
  parser,
  opt_str = c("--hg_sets"), type = "character", action = "store",
  help = "COMMA separated homology group sets"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-c", "--config"), type = "character", action = "store",
  help = "project config YAML file"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-o", "--out"), type = "character", action = "store",
  default = "HG_pav.pdf", help = "Output file name"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$hg_sets), is.null(opts$config), is.null(opts$tree))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# ## for test
# opts$config <- "project_config.yaml"
# opts$tree <- "kmer_upgma"
# opts$out <- "assay_HG_pav"
# opts$hg_sets <- "LZI,TIR,assay_FN"
################################################################################

setIds <- stringr::str_split(opts$hg_sets, pattern = ",") %>%
  unlist()

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = opts$config)),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

# add any homology group set file created to this vector
# `setIds` are searched in these files to identify respective homology groups
# for plotting
hgSetFiles <- c(
  confs$analysis$homology_groups$files$hg_sets,
  confs$analysis$association$files$pheno_specific_groups,
  confs$analysis$prophages$preprocessing$files$raw_prophage_hg
)

################################################################################


rawTree <- import_tree(
  file = confs$analysis$phylogeny[[opts$tree]]$files$tree, phylo = TRUE
)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny[[opts$tree]]$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) %>%
  dplyr::mutate(
    SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
  )

################################################################################

# prepare homology groups data for plotting
hgSets <- NULL

for (st in hgSetFiles) {
  df <- suppressMessages(readr::read_tsv(st)) %>%
    dplyr::mutate(
      hgs = stringr::str_split(hgs, ";"),
      setSource = st
    ) %>%
    dplyr::select(id, hgs, setSource)


  hgSets <- dplyr::bind_rows(hgSets, df)
}



if (any(duplicated(hgSets$id))) {
  stop(
    "duplicate homology group set identifiers found in the data:",
    names(which(table(hgSets$id) > 1))
  )
}

# read prophage HGs stored locally
hgL <- dplyr::left_join(
  x = tibble::tibble(id = setIds), y = hgSets, by = "id"
) %>%
  purrr::transpose() %>%
  purrr::set_names(nm = purrr::map(., "id"))


# species key heatmap
htList <- species_key_heatmap(
  genomes = rawTree$tip.label, speciesInfo = sampleInfo,
  use_raster = TRUE, raster_quality = 3
)

htList@heatmap_param$width <- unit(12, "cm")

# prepare homology group PAV matrix from pan.db
hgPavMat <- homology_groups_mat(
  pandb = panOrgDb, type = "pav",
  groups = purrr::map(hgL, "hgs") %>% unlist(use.names = F) %>% unique()
)

hgPavMat <- hgPavMat[rawTree$tip.label, ]

for (id in names(hgL)) {

  hgMat <- hgPavMat[, hgL[[id]]$hgs, drop = FALSE]

  ht <- homology_group_heatmap(
    mat = hgMat, phy = rawTree,
    # width = unit(10, "cm"),
    name = id, column_title = id,
    use_raster = TRUE, raster_quality = 3
  )

  ht@column_dend_param$cluster <- FALSE
  ht@column_names_param$show <- FALSE

  htList <- htList + ht
}


pdf(file = opts$out, width = 12, height = 8)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = hgL[[1]]$id,
  row_dend_side = "left",
  merge_legends = TRUE,
  heatmap_legend_side = "bottom"
)
dev.off()
