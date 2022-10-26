suppressPackageStartupMessages(library(tidyverse))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## prepare input data for pangenome construction using PanTools

rm(list = ls())

#####################################################################

file_metadata <- here::here("data/reference_data", "sample_metadata.tsv")
file_duplicateGenomes <- here::here("analysis", "01_multiqc", "duplicate_genomes.tab")
file_excludeGenomes <- here::here("analysis", "01_multiqc", "exclude_genomes.txt")

path_genomes <- here::here("data", "prokka_annotation")

cutoff_buscog <- 99

#####################################################################

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
    )
  ) %>% 
  dplyr::relocate(filtered, .after = AssemblyName)

filteredMeta <- dplyr::filter(metadata, filtered == "PASS") %>% 
  dplyr::mutate(
    genomeId = 1:n(),
    fasta = paste(path_genomes, "/", sampleId, "/", sampleId, ".fna", sep = ""),
    gff3 = paste(path_genomes, "/", sampleId, "/", sampleId, ".gff", sep = "")
  )


## FASTA file paths
dplyr::select(filteredMeta, fasta) %>% 
  readr::write_tsv(
    file = here::here("analysis", "04_pangenome", "genomes_fa.list"),
    col_names = FALSE
  )

## GFF3 file paths
dplyr::select(filteredMeta, genomeId, gff3) %>% 
  readr::write_tsv(
    file = here::here("analysis", "04_pangenome", "genomes_gff3.list"),
    col_names = FALSE
  )

## metadata file
dplyr::select(
  filteredMeta, genomeId, sampleId, AssemblyAccession,	AssemblyName,	SpeciesName,
  taxonomy_check_status, strain, virulence, virulence_pcr, geo_loc_country, host,
  isolation_source,	collected_by, env_broad_scale, type_material, virulence, virulence_pcr
  ) %>% 
  readr::write_csv(
    file = here::here("analysis", "04_pangenome", "genomes_metadata.csv"),
    col_names = FALSE
  )




