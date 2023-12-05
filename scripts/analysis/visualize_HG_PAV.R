#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

## plot PAV for the homology groups specific for a particular phenotype

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/GO_enrichment/enrichment_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/association_analysis.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

phenotype <- "assay_FN"
treeMethod <- "ani_upgma"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- file.path(confs$analysis$association$dir, 'pyseer_clades', phenotype)
outPrefix <- file.path(outDir, phenotype)
file_result <- file.path(outDir, 'pav.virulent_Pbrasiliense.asso.no_dist.significant.txt')
orgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## genomes that are used for phenotype association analysis
associatedGenomes <- get_phenotype_association_genomes(
  phenotype = phenotype,
  confFile = panConf$analysis_confs$files$clade_association
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>%  dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

## import homology group data
g2hg <- suppressMessages(
  readr::read_tsv(confs$analysis$homology_groups$files$groups)
) %>% 
  dplyr::mutate(Genome = as.character(Genome))

assoRes <- suppressMessages(readr::read_tsv(file_result)) %>% 
  dplyr::filter(is.na(notes)) %>%
  dplyr::mutate(
    hg = stringr::str_replace(variant, "hg_", ""),
    betaSign = sign(beta),
    interceptSign = sign(intercept)
  )


################################################################################
hgMat <- dplyr::select(g2hg, Genome, assoRes$hg) %>% 
  tibble::column_to_rownames(var = "Genome") %>% 
  as.matrix()

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo,
  hgAn = assoRes,
  width = c(10,24),
  markGenomes = associatedGenomes
)

png(filename = paste(outPrefix, ".pheno_hg_association.png", sep = ""), width = 6000, height = 3000, res = 350)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
)
dev.off()

################################################################################

# topGoDf <- topGO_enrichment(genes = assoRes$hg, orgdb = orgDb)

## GO enrichment 
grpGo <- dplyr::group_by(assoRes, betaSign) %>% 
  dplyr::group_modify(
    .f = ~topGO_enrichment(genes = .x$hg, orgdb = orgDb)
  ) %>% 
  dplyr::ungroup()

readr::write_tsv(
  grpGo,
  file = paste(outPrefix, ".topGO.tab", sep = "")
)

(pt_go <- enrichment_bar(df = grpGo, title = "GO", colorCol = "betaSign"))

ggsave(
  plot = pt_go, filename = paste(outPrefix, ".topGO.png", sep = ""),
  width = 8, height = 5
)

