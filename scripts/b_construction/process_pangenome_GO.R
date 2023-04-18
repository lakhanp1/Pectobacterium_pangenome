suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(AnnotationForge))

## compare homology groups counts between previous and current Pecto pangenomes
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

outDir <- panConf$dir
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata) %>% 
  dplyr::select(Genome, SpeciesName)

panDf <- suppressMessages(
  readr::read_tsv(file = panConf$files$go_data, na = "None")
) %>% 
  dplyr::rename(GID = hg_id) %>% 
  dplyr::mutate(GID = as.character(GID))

hgDf <- dplyr::select(panDf, GID , mRNA_id, genome, chr) %>% 
  dplyr::distinct() %>% 
  tidyr::unite(genome, chr, mRNA_id, col = "mRNA_key", sep = "_", remove = FALSE) %>% 
  dplyr::select(GID, everything())

hgGoDf <- dplyr::select(panDf, GID, GO = go_id) %>% 
  dplyr::filter(!is.na(GO)) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(EVIDENCE = "IEA")

hgMeta <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups
  )
) %>% 
  dplyr::select(GID = "Homology group id", class) %>% 
  dplyr::mutate(GID = as.character(GID))

#sed -n '1,/##FASTA/p' interproscan/GCF_024506455.1_ASM2450645v1.interProScan.gff3 | grep -v '^#' | cut -f1 | sort | uniq | wc -l
#sed -n '1,/##FASTA/p' interproscan/GCF_024506455.1_ASM2450645v1.interProScan.gff3 | grep -v '^#' | grep 'GO:' | cut -f1 | sort | uniq | wc -l
################################################################################

## prepare pangenome.db


AnnotationForge::makeOrgPackage(
  hgDf = hgDf,
  hgMeta = hgMeta,
  go = hgGoDf,
  version = "1.0.0",
  maintainer = "Lakhansing Pardeshi <lakhansing.pardeshi@wur.nl>",
  author = "Lakhansing Pardeshi",
  outputDir = "./data",
  tax_id = "122277",
  genus = "P",
  species = "ectobacterium.spp.pan",
  goTable = "go",
  verbose = TRUE
)

devtools::install(
  pkg = "./data/org.Pectobacterium.spp.pan.eg.db",
  upgrade = "never"
)

orgDb <- org.Pectobacterium.spp.pan.eg.db
keys(orgDb)
