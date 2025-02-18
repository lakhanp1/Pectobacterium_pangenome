suppressPackageStartupMessages(library(tidyverse))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## prepare input data for pangenome construction using PanTools

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

panConf <- confs$data$pangenomes$pectobacterium.v2
testPanConf <- confs$data$pangenomes$pectobacterium.ts

cols_metadata <- c(
  "sampleName", "AssemblyAccession", "AssemblyName", "SpeciesName", "taxonomy_check_status",
  "BioSampleAccn", "BioprojectAccn",
  "strain", "virulence", "virulence_pcr", "geo_loc_country", "continent",
  "host", "isolation_source", "collection_year", "collected_by", "env_broad_scale",
  "type_material", "virulence", "virulence_pcr",
  "source", "type_material", "representative_status", "sample_type",
  "length", "N50", "L50", "n_contigs"
)

#####################################################################
!dir.exists(panConf$path) && dir.create(path = panConf$path, recursive = TRUE)
!dir.exists(testPanConf$path) && dir.create(path = testPanConf$path, recursive = TRUE)

if (!dir.exists(panConf$path)) {
  dir.create(path = panConf$path, recursive = TRUE)
}

metadata <- suppressMessages(
  readr::read_tsv(file = confs$data$reference_data$files$metadata)
)


filteredMeta <- dplyr::filter(
  metadata,
  filtered == "PASS",
  source %in% panConf$include_source
) %>%
  dplyr::mutate(
    Genome = 1:n(),
    genomeId = paste("g_", Genome, sep = ""),
    fasta = paste(
      confs$data$prokka$path, "/", sampleId, "/", sampleId, ".fna",
      sep = ""
    ),
    gff3 = paste(
      confs$data$prokka$path, "/", sampleId, "/", sampleId, ".gff3",
      sep = ""
    ),
    interpro = paste(
      confs$data$interproscan$path, "/", sampleId, ".interProScan.gff3",
      sep = ""
    ),
    cog = paste(
      confs$data$cog$path, "/", sampleId, ".emapper.annotations",
      sep = ""
    ),
    ## replace "," with ";"
    dplyr::across(
      .cols = everything(),
      .fns = ~ stringr::str_replace_all(
        string = .x, pattern = "[,;]", replacement = " and"
      )
    )
  )


#####################################################################
## FASTA file paths
dplyr::select(filteredMeta, fasta) %>%
  readr::write_tsv(
    file = panConf$files$genomes,
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(filteredMeta, Genome, gff3) %>%
  readr::write_tsv(file = panConf$files$gff, col_names = FALSE)

## functional annotation file paths
dplyr::select(filteredMeta, Genome, interpro) %>%
  readr::write_delim(file = panConf$files$annotations, col_names = FALSE)

## COG annotation file path
dplyr::select(filteredMeta, Genome, cog) %>%
  readr::write_delim(file = panConf$files$cog, col_names = FALSE)

## metadata file
dplyr::select(filteredMeta, Genome, genomeId, sampleId, !!!cols_metadata) %>%
  readr::write_csv(file = panConf$files$metadata, col_names = TRUE)

## input genome lock file: this file should never change
readr::write_tsv(
  x = dplyr::select(filteredMeta, Genome, sampleId),
  file = panConf$files$input_lock
)

## write metadata to excel
wb <- openxlsx::createWorkbook()
currentSheet <- "metadata"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, withFilter = TRUE, keepNA = TRUE, na.string = "NA",
  x = dplyr::select(
    filteredMeta,
    Genome = genomeId, sampleId, !!!cols_metadata
  )
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = file.path(panConf$path, "pangenome_metadata.xlsx"), overwrite = TRUE
)


#####################################################################
## write small subset for testing pangenome pipeline
testSet <- dplyr::filter(filteredMeta, !is.na(type_material)) %>%
  dplyr::arrange(sampleId) %>%
  dplyr::mutate(
    Genome = 1:n(),
    genomeId = paste("g_", Genome, sep = "")
  )

## FASTA file paths
dplyr::select(testSet, fasta) %>%
  readr::write_tsv(file = testPanConf$files$genomes, col_names = FALSE)

## GFF3 file paths
dplyr::select(testSet, Genome, gff3) %>%
  readr::write_tsv(file = testPanConf$files$gff, col_names = FALSE)

## functional annotation file paths
dplyr::select(testSet, Genome, interpro) %>%
  readr::write_delim(file = testPanConf$files$annotations, col_names = FALSE)

## metadata file
dplyr::select(testSet, Genome, sampleId, !!!cols_metadata) %>%
  readr::write_csv(file = testPanConf$files$metadata, col_names = TRUE)
