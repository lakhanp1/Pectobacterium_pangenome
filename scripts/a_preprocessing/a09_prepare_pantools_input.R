suppressPackageStartupMessages(library(tidyverse))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## prepare input data for pangenome construction using PanTools

rm(list = ls())
set.seed(124)
#####################################################################

file_metadata <- here::here("data/reference_data", "sample_metadata.tsv")
file_duplicateGenomes <- here::here("analysis", "01_multiqc", "duplicate_genomes.tab")
file_excludeGenomes <- here::here("analysis", "01_multiqc", "exclude_genomes.txt")

pangenomeName <- "pectobacterium.v2"
testPangenome <- "pectobacterium.10g"
path_genomes <- "./data/prokka_annotation"
path_out <- here::here("data", "pangenomes", pangenomeName)
test_out <- here::here("data", "pangenomes", testPangenome)

cutoff_buscog <- 99

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

metadata <- suppressMessages(readr::read_tsv(file = file_metadata))
duplicateGenomes <- suppressMessages(readr::read_tsv(file = file_duplicateGenomes)) %>% 
  dplyr::filter(identical == TRUE)
excludeGenomes <- suppressMessages(
  readr::read_tsv(file = file_excludeGenomes, col_names = "sampleId")
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
    fasta = paste(path_genomes, "/", sampleId, "/", sampleId, ".fna", sep = ""),
    gff3 = paste(path_genomes, "/", sampleId, "/", sampleId, ".gff", sep = ""),
    interpro = paste(
      "./data/interproscan/", sampleId, ".interProScan.gff3", sep = ""
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
    file = file.path(path_out, "genomes_fa.list"),
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(filteredMeta, genomeId, gff3) %>% 
  readr::write_tsv(
    file = file.path(path_out, "genomes_gff3.list"),
    col_names = FALSE
  )

## functional annotation file paths
dplyr::select(filteredMeta, genomeId, interpro) %>% 
  readr::write_delim(
    file = file.path(path_out, "functional_annotations.txt"),
    col_names = FALSE
  )

## metadata file
dplyr::select(
  filteredMeta, Genome=genomeId, id=sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = file.path(path_out, "genomes_metadata.csv"),
    col_names = TRUE
  )

## input genome lock file: this file should never change
readr::write_tsv(
  x = dplyr::select(filteredMeta, Genome = genomeId, id = sampleId),
  file = file.path(path_out, "input_genomes.tab")
)

#####################################################################
## write small subset for testing pangenome pipeline
testSet <- dplyr::filter(filteredMeta, SpeciesName == "Pectobacterium brasiliense") %>% 
  dplyr::slice_sample(n = 10) %>% 
  dplyr::arrange(sampleId) %>% 
  dplyr::mutate(genomeId = 1:n())

## FASTA file paths
dplyr::select(testSet, fasta) %>% 
  readr::write_tsv(
    file = file.path(test_out, "genomes_fa.list"),
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(testSet, genomeId, gff3) %>% 
  readr::write_tsv(
    file = file.path(test_out, "genomes_gff3.list"),
    col_names = FALSE
  )

## functional annotation file paths
dplyr::select(testSet, genomeId, interpro) %>% 
  readr::write_delim(
    file = file.path(test_out, "functional_annotations.txt"),
    col_names = FALSE
  )

## metadata file
dplyr::select(
  testSet, Genome=genomeId, id=sampleId, !!!cols_metadata
) %>% 
  readr::write_csv(
    file = file.path(test_out, "genomes_metadata.csv"),
    col_names = TRUE
  )

