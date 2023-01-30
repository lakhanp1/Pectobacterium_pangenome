#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))


## plot default tree with node lables information for downstream analysis
## inspection

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
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

analysisName <- "default_tree"
outDir <- confs$analysis$phylogeny$dir
outPrefix <- file.path(outDir, analysisName)

################################################################################

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## read tree
rawTree <- ape::read.tree(confs$analysis$phylogeny$files$default_tree)

## set negative length edges => 0
rawTree$edge.length[rawTree$edge.length < 0] <- 0

rootedTr <- ape::root(phy = rawTree, outgroup = sampleInfoList[[outGroup]]$Genome) %>% 
  ape::ladderize()

## add data to tree
treeTbl <- as_tibble(rootedTr) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata()

################################################################################

pt_tree <- ggtree::ggtree(tr = treeTbl) +
  labs(title = "kmer distance NJ tree")

## mark outgroup
pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleName")

## mark species of interest
pt_tree3 <- pt_tree2 +
  ggtree::geom_nodelab(
    mapping = aes(label = node),
    node = "internal", size = 3, hjust = 1.3
  ) +
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


ggsave(
  plot = pt_tree4, width = 10, height = 20, scale = 2,
  filename = paste(outPrefix, "_tree.pdf", sep = "")
)

# ggtree::viewClade(tree_view = pt_tree4, node = 626)
# ggtree::scaleClade(tree_view = pt_tree4, node = 626, scale = 8)

subTree2 <- ape::extract.clade(phy = rawTree, node = 626)
plot(ladderize(subTree2))

## group genomes from clade of interest into a phenotype category
ph1 <- prepare_comparison_data(tr = treeTbl, node = 563, name = "Pbrasiliense_clade")
ph2 <- prepare_comparison_data(
  tr = treeTbl, node = 626, zoom = 563, name = "virulent_Pbrasiliense"
)
ph3 <- prepare_comparison_data(
  tr = treeTbl, node = 661, zoom = 563, name = "assay_FN", category = "assay_FN"
)

## write updated phenotype file with clade of interest as phenotypes
dplyr::select(sampleInfo, Genome, everything(), -nodeLabs) %>% 
  dplyr::left_join(y = ph1$pheno, by = "Genome") %>% 
  dplyr::left_join(y = ph2$pheno, by = "Genome") %>% 
  dplyr::left_join(y = ph3$pheno, by = "Genome") %>% 
  dplyr::mutate(
    dplyr::across(
      .fns = ~stringr::str_replace_all(.x, ";", " and")
    )
  ) %>% 
  readr::write_csv(
    file = paste(outDir, "/clade_comarison_phenotypes.csv", sep = "")
  )

## `pantools gene_classification` configuration
list(ph1, ph2, ph3) %>% 
  purrr::map_dfr(
    .f = function(x){
      return(list(name = x$name, compare = x$compare, against = x$against))
    }
  ) %>% 
  dplyr::mutate(
    phenotypeArg = dplyr::if_else(
      condition = is.na(against) | against == "",
      true = stringr::str_c("--phenotype=", name, sep = ""), 
      false = stringr::str_c("--phenotype=", name, " --include=", against, sep = "")
    )
  ) %>% 
  dplyr::select(name, phenotypeArg, compare, against) %>% 
  readr::write_tsv(
    file = confs$data$analysis_confs$files$phenotype_association,
    col_names = FALSE, na = ""
  )
