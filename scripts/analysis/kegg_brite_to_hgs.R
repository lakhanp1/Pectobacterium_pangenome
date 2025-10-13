#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
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
panOrgDb <- org.Pectobacterium.spp.pan.eg.db

# brightId <- "01003" # Glycosyltransferases
brightId <- "01005" # Lipopolysaccharide Biosynthesis Proteins
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) |>
  dplyr::mutate(
    GenBank = stringr::str_remove(AssemblyAccession, "\\.\\d+$")
  )

################################################################################

# get required genomes from KEGG: https://rest.kegg.jp/list/genome
org <- as.data.frame(keggList("organism")) |>
  dplyr::filter(stringr::str_detect(species, "^Pectobacterium")) |>
  dplyr::mutate(
    genomePage = paste("https://rest.kegg.jp/get/gn:", organism, sep = ""),
    genesPage = paste("https://rest.kegg.jp/list/", organism, sep = "")
  )


i <- 1

orgCode <- org$organism[i]

# get Lipopolysaccharide Biosynthesis Proteins from BRIGHT hierarchy
# https://rest.kegg.jp/list/brite
# https://rest.kegg.jp/get/br:ko01003/json
# https://rest.kegg.jp/get/br:pbra01005/json
# keggGet("br:ko01005")
bright <- keggGet(paste("br:", orgCode, brightId, sep = ""), "json")

# parse BRITE proteins JSON output
briteProteins <- jsonlite::fromJSON(txt = bright)$children |>
  dplyr::mutate(
    org = .env$orgCode,
    children = purrr::map(children, "name")
  ) |>
  tidyr::unnest(cols = children) |>
  tidyr::separate_wider_regex(
    cols = children,
    pattern = c(gene = "\\w+", " ", description = ".*", "\\t", ko = "\\w+", " ",
                ko_symbol = ".*?", "; ", ko_name = ".*")
  )

# https://rest.kegg.jp/get/gn:pvz
genome <- keggGet(paste("gn:", orgCode, sep = ""))[[1]]

# https://rest.kegg.jp/list/pvz
genes <- suppressMessages(readr::read_tsv(
  file = org$genesPage[i],
  col_names = c("gene", "type", "position", "description")
)) |>
  tidyr::separate_wider_delim(cols = gene, delim = ":", names = c("org", "gene")) |>
  dplyr::select(-description) |>
  dplyr::mutate(
    strand = stringr::str_detect(position, "^complement\\("),
    strand = if_else(strand, "-", "+"),
    start = stringr::str_extract(position, pattern = "(\\d+)\\.\\.", group = 1),
    end = stringr::str_extract(position, pattern = "\\.\\.(\\d+)", group = 1)
  )

briteProteins <- dplyr::left_join(
  briteProteins, genes, by = c("org", "gene")
)
