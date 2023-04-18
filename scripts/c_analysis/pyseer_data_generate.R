#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))

## generate input files for pyseer association analysis of virulence in Pbr isolates

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

analysisName <- "pbr_virulence"
outDir <- file.path(confs$analysis$association$dir, analysisName)
outPrefix <- file.path(outDir, analysisName)
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata)
rawTree <- ape::read.tree(file = panConf$db$core_phylogeny$default$files$tree_ml)
rawTree <- ape::updateLabel(rawTree, old = sampleInfo$Genome, new = sampleInfo$sampleId, exact = TRUE)

## input for unitig-counter 
sampleInfo %<>% dplyr::filter(SpeciesName == "P. brasiliense") %>% 
  dplyr::mutate(
    fasta = paste(
      confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""
    )
  )

if(!dir.exists(outDir)){
  dir.create(outDir)
}

readr::write_lines(
  sampleInfo$fasta,
  file = paste(outDir, "/", analysisName, ".fa.list", sep = "")
)

## homology groups PAV matrix


## cophenetic distance of phylogenetic tree
distMat <- ape::cophenetic.phylo(rawTree)

distMat <- distMat[sampleInfo$sampleId, sampleInfo$sampleId]
write.table(
  x = distMat, col.names = NA, row.names = TRUE, sep = "\t", quote = F,
  file = paste(outDir, "/phylogeny_dist.tab", sep = "")
)


## save the phenotype
dplyr::select(sampleInfo, samples = sampleId, Genome, virulence) %>% 
  dplyr::mutate(
    virulence = dplyr::if_else(
      condition = virulence == "virulent", true = 1, false = 0
    )
  ) %>% 
  readr::write_tsv(file = paste(outDir, "/phenotypes.tab", sep = ""))

## save the subtree for similarity calculation
pbrTree <- ape::keep.tip(rawTree, tip = sampleInfo$sampleId)
ape::write.tree(
  phy = pbrTree,
  file = paste(outDir, "/phylogeny_tree.newick", sep = "")
)

