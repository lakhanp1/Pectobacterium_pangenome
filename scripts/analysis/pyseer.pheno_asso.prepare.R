#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(matrixStats))

## generate input files for pyseer association analysis of virulence in Pbr isolates

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

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

analysisName <- "pbr_virulence2"
genomeSetName <- "pbr"
outDir <- file.path(confs$analysis$association$path, analysisName)
outPrefix <- file.path(outDir, analysisName)
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)
rawTree <- ape::read.tree(file = panConf$db$core_phylogeny$default$files$tree_ml)
rawTree <- ape::updateLabel(rawTree, old = sampleInfo$Genome, new = sampleInfo$sampleId, exact = TRUE)

## input for unitig-counter
sampleInfo %<>% dplyr::filter(SpeciesName == "P. brasiliense") %>%
  dplyr::filter(!is.na(virulence)) %>%
  dplyr::mutate(
    fasta = paste(
      confs$data$prokka$path, "/", sampleId, "/", sampleId, ".fna",
      sep = ""
    )
  )

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

if (!dir.exists(outDir)) {
  dir.create(outDir)
}
################################################################################
readr::write_lines(
  sampleInfo$fasta,
  file = paste(outDir, "/", genomeSetName, ".fa.list", sep = "")
)

## homology groups
hgs <- homology_groups_extract(
  file = panConf$db$gene_classification$GC.100.0$files$groups,
  genomes = sampleInfo$Genome,
  groups <- "accessory", pav = TRUE
)

dplyr::select(hgs, -class, Gene = hg, !!!purrr::map_chr(sampleInfoList, "Genome")) %>%
  dplyr::mutate(Gene = paste("hg_", Gene, sep = "")) %>%
  readr::write_tsv(
    file = paste(outDir, "/", genomeSetName, ".accessory_PAV.tab", sep = "")
  )


## cophenetic distance of phylogenetic tree
distMat <- ape::cophenetic.phylo(rawTree)

distMat <- distMat[sampleInfo$sampleId, sampleInfo$sampleId]
write.table(
  x = distMat, col.names = NA, row.names = TRUE, sep = "\t", quote = F,
  file = paste(outDir, "/phylogeny_dist.tab", sep = "")
)

mds <- cmdscale(d = distMat, k = 5, eig = TRUE)

ggplot(data = tibble(eig = mds$eig, n = 1:length(mds$eig))) +
  geom_bar(mapping = aes(x = n, y = eig), stat = "identity")

tibble::as_tibble(x = mds$points, rownames = "sampleId") %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, Genome, virulence),
    by = "sampleId"
  ) %>%
  tidyr::replace_na(replace = list(virulence = "NA")) %>%
  # view() %>%
  ggplot(mapping = aes(x = V1, y = V2, fill = virulence, label = virulence)) +
  geom_point(size = 3, alpha = 0.8, color = "black", shape = 22) +
  ggrepel::geom_text_repel(size = 2)

## save the phenotype
dplyr::select(sampleInfo, samples = sampleId, Genome, virulence) %>%
  dplyr::mutate(
    virulence = dplyr::if_else(
      condition = virulence == "virulent", true = 1, false = 0
    )
  ) %>%
  readr::write_tsv(file = paste(outDir, "/", genomeSetName, "_phenotypes.tab", sep = ""))

## save the subtree for similarity calculation
pbrTree <- ape::keep.tip(rawTree, tip = sampleInfo$sampleId)
ape::write.tree(
  phy = pbrTree,
  file = paste(outDir, "/phylogeny_tree.newick", sep = "")
)
