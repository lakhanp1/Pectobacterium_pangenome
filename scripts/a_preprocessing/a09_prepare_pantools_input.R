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

pangenomeName <- confs$data$pangenomes$pectobacterium.v2$name
testPangenome <- confs$data$pangenomes$pectobacterium.10g$name
path_out <- confs$data$pangenomes[[pangenomeName]]$dir
test_out <- confs$data$pangenomes[[testPangenome]]$dir

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

filteredMeta <- dplyr::filter(
  metadata,
  filtered == "PASS",
  source %in% confs$data$include_source
) %>% 
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
        string = .x, pattern = "[,;]", replacement = " and"
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
  filteredMeta, Genome=genomeId, sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = confs$data$pangenomes[[pangenomeName]]$files$metadata,
    col_names = TRUE
  )

## input genome lock file: this file should never change
readr::write_tsv(
  x = dplyr::select(filteredMeta, Genome = genomeId, sampleId),
  file = confs$data$pangenomes[[pangenomeName]]$files$input_lock
)

## write metadata to excel
wb <- openxlsx::createWorkbook()
currentSheet <- "metadata"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, withFilter = TRUE, keepNA = TRUE, na.string = "NA",
  x = dplyr::select(
    filteredMeta, Genome=genomeId, sampleId, !!!cols_metadata
  )
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = file.path(path_out, "pangenome_metadata.xlsx"), overwrite = TRUE
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
  testSet, Genome=genomeId, sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = confs$data$pangenomes[[testPangenome]]$files$metadata,
    col_names = TRUE
  )

