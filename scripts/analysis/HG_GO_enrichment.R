#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# GO enrichment for homology groups

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/GO_enrichment/enrichment_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

analysisName <- "homology_groups"
outDir <- confs$analysis$homology_groups$dir
outPrefix <- file.path(outDir, analysisName)
panOrgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

spOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$core_snp_ml$files$species_order)
)

hgTable <- AnnotationDbi::select(
  x = panOrgDb, keys = keys(panOrgDb),
  columns = c("GID", "genomeId", "class")
) %>%
  dplyr::rename(hgId = GID) %>%
  dplyr::count(hgId, genomeId, name = "nGenes")

hgMeta <- AnnotationDbi::select(
  x = panOrgDb, keys = keys(panOrgDb),
  columns = c("GID", "class")
) %>%
  dplyr::rename(hgId = GID) %>%
  dplyr::mutate(
    class = dplyr::if_else(
      condition = class == "core & single copy orthologous",
      true = "core", false = class
    )
  )

## binary matrix for homology_group PAV
hgPavMat <- homology_groups_mat(pandb = panOrgDb, type = "pav")

################################################################################
## GO enrichment for homology groups

spNames <- dplyr::count(sampleInfo, SpeciesName) %>%
  # dplyr::filter(n >= 20) %>%
  dplyr::pull(SpeciesName)

sppGrpStats <- NULL
sppGrpGo <- NULL


## get species wise core, accessory, unique group stats and GO
for (sp in spNames) {
  spGenomes <- dplyr::filter(sampleInfo, SpeciesName == .env$sp) %>%
    dplyr::pull(genomeId)

  cat(sp, ": ", length(spGenomes), "\n")

  hgSum <- matrixStats::colSums2(
    x = hgPavMat, useNames = T,
    rows = which(rownames(hgPavMat) %in% spGenomes)
  ) %>%
    tibble::enframe(name = "hgId", value = "nGenomes") %>%
    dplyr::filter(nGenomes != 0) %>%
    dplyr::mutate(
      subpan_class = dplyr::case_when(
        nGenomes == !!length(spGenomes) ~ "core",
        nGenomes == 1 ~ "unique",
        nGenomes < !!length(spGenomes) & nGenomes > 1 ~ "accessory"
      )
    ) %>%
    dplyr::left_join(
      y = dplyr::rename(hgMeta, pangenome_class = class),
      by = "hgId"
    )

  ## GO enrichment
  sppGrpGo <- dplyr::group_by(hgSum, subpan_class) %>%
    dplyr::group_modify(
      .f = ~ topGO_enrichment(genes = .x$hgId, orgdb = panOrgDb)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SpeciesName = .env$sp) %>%
    dplyr::bind_rows(sppGrpGo)
}


panGo <- hgMeta  %>%
  dplyr::group_by(class) %>%
  dplyr::group_modify(
    .f = ~ topGO_enrichment(genes = .x$hgId, orgdb = panOrgDb)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SpeciesName = "pangenome")

sppGrpGo <- dplyr::bind_rows(panGo, sppGrpGo) %>%
  dplyr::select(SpeciesName, everything())

readr::write_tsv(
  x = sppGrpGo, file = confs$analysis$homology_groups$files$spp_group_go
)
