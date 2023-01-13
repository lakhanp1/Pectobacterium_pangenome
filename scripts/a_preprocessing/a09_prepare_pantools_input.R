suppressPackageStartupMessages(library(tidyverse))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## prepare input data for pangenome construction using PanTools

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

# 
# file_metadata <- here::here("data/reference_data", "sample_metadata.tsv")
# file_duplicateGenomes <- here::here("analysis", "01_QC", "duplicate_genomes.tab")
# file_excludeGenomes <- here::here("analysis", "01_QC", "exclude_genomes.txt")

pangenomeName <- confs$data$pangenomes$pectobacterium.v2$name
testPangenome <- confs$data$pangenomes$pectobacterium.10g$name
path_out <- confs$data$pangenomes[[pangenomeName]]$dir
test_out <- confs$data$pangenomes[[testPangenome]]$dir

cutoff_buscog <- confs$parameters$cutoff_busco

cols_metadata <- c(
  "sampleName", "AssemblyAccession",	"AssemblyName", "SpeciesName", "taxonomy_check_status", 
  "strain", "virulence", "virulence_pcr", "geo_loc_country", "host", "isolation_source",
  "collection_year", "collected_by", "env_broad_scale", 
  "type_material", "virulence", "virulence_pcr",
  "source", "type_material", "representative_status", "sample_type")

#####################################################################
! dir.exists(path_out) && dir.create(path = path_out, recursive = TRUE)
! dir.exists(test_out) && dir.create(path = test_out, recursive = TRUE)

if(! dir.exists(path_out)){
  dir.create(path = path_out, recursive = TRUE)
}

metadata <- suppressMessages(
  readr::read_tsv(file = confs$data$reference_data$files$metadata)
)
duplicateGenomes <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$duplicate_genomes)
) %>% 
  dplyr::filter(identical == TRUE)
excludeGenomes <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$exclude_genomes, col_names = "sampleId")
)


metadata %<>% 
  dplyr::mutate(
    filtered = dplyr::case_when(
      !is.na(ExclFromRefSeq) ~ "ExclFromRefSeq",
      !is.na(Anomalous) ~ "Anomalous",
      !is.na(replaced) ~ "replaced",
      sampleId %in% c(duplicateGenomes$genome2, excludeGenomes$sampleId) ~ "duplicate",
      buscog.complete < cutoff_buscog ~ paste("BUSCO <", cutoff_buscog),
      TRUE ~ "PASS"
    ),
    type_material = stringr::str_replace(
      string = tolower(type_material), pattern = "(type (strain|material)).*", 
      replacement = "type strain"
    )
  ) %>% 
  dplyr::relocate(filtered, .after = AssemblyName)

filteredMeta <- dplyr::filter(metadata, filtered == "PASS") %>% 
  dplyr::mutate(
    genomeId = 1:n(),
    fasta = paste(
      confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""
    ),
    gff3 = paste(
      confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".gff", sep = ""
    ),
    interpro = paste(
      confs$data$interproscan$dir, "/", sampleId, ".interProScan.gff3", sep = ""
    ),
    ## replace "," with ";"
    dplyr::across(
      .cols = everything(),
      .fns = ~stringr::str_replace_all(
        string = .x, pattern = ",", replacement = ";"
      )
    )
  )

#####################################################################
## FASTA file paths
dplyr::select(filteredMeta, fasta) %>% 
  readr::write_tsv(
    file = confs$data$pangenomes[[pangenomeName]]$files$genomes,
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(filteredMeta, genomeId, gff3) %>% 
  readr::write_tsv(
    file = confs$data$pangenomes[[pangenomeName]]$files$gff,
    col_names = FALSE
  )

## functional annotation file paths
dplyr::select(filteredMeta, genomeId, interpro) %>% 
  readr::write_delim(
    file = confs$data$pangenomes[[pangenomeName]]$files$annotations,
    col_names = FALSE
  )

## metadata file
dplyr::select(
  filteredMeta, Genome=genomeId, id=sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = confs$data$pangenomes[[pangenomeName]]$files$metadata,
    col_names = TRUE
  )

## input genome lock file: this file should never change
readr::write_tsv(
  x = dplyr::select(filteredMeta, Genome = genomeId, id = sampleId),
  file = confs$data$pangenomes[[pangenomeName]]$files$input_lock
)

#####################################################################
## write small subset for testing pangenome pipeline
testSet <- dplyr::filter(filteredMeta, !is.na(type_material)) %>% 
  dplyr::arrange(sampleId) %>% 
  dplyr::mutate(genomeId = 1:n())

## FASTA file paths
dplyr::select(testSet, fasta) %>% 
  readr::write_tsv(
    file = confs$data$pangenomes[[testPangenome]]$files$genomes,
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(testSet, genomeId, gff3) %>% 
  readr::write_tsv(
    file = confs$data$pangenomes[[testPangenome]]$files$gff,
    col_names = FALSE
  )

## functional annotation file paths
dplyr::select(testSet, genomeId, interpro) %>% 
  readr::write_delim(
    file = confs$data$pangenomes[[testPangenome]]$files$annotations,
    col_names = FALSE
  )

## metadata file
dplyr::select(
  testSet, Genome=genomeId, id=sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = confs$data$pangenomes[[testPangenome]]$files$metadata,
    col_names = TRUE
  )

