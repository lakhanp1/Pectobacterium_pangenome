suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))

## check the taxonomy for bacterial genomes
## generate a table for manual inspection

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- confs$analysis$qc$dir

################################################################################

genomeMetadata <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$prebuild_metadata)
)

ncbiAni <- suppressMessages(
  readr::read_tsv(file = confs$data$other$files$ani_report, na = "na")
)

typeStrains <- dplyr::filter(genomeMetadata, !is.na(type_material)) %>% 
  dplyr::select(type_organism_name = SpeciesName, typeStrainId = sampleId)

aniDf <- suppressMessages(readr::read_tsv(
  file = confs$analysis$ANI$files$fastani_out,
  col_names = c("id1", "id2", "ani", "mapped", "total")
)) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = c(id1, id2),
      .fns = ~stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1"),
    )
  )

typeStrainAni <- dplyr::left_join(
  typeStrains, aniDf, by = c("typeStrainId" = "id2")
)

################################################################################
## taxonomy check failed/inconclusive information

## NCBI data
ncbiTaxCheck <- dplyr::filter(genomeMetadata, source == "NCBI") %>%
  dplyr::select(
    sampleId, source, AssemblyAccession, synonym_GB, taxonomy_check_status
  ) %>%
  dplyr::filter(
    taxonomy_check_status != "OK"
  ) %>%
  dplyr::left_join(
    y = dplyr::select(
      ncbiAni, genbank_accession, species_name, organism_name,
      declared_type_organism_name, declared_type_ANI,
      declared_type_qcoverage, declared_type_scoverage, best_match_species_name,
      best_match_type_ANI, best_match_type_qcoverage, best_match_type_scoverage,
      best_match_status, comment
    ),
    by = c("synonym_GB" = "genbank_accession")
  ) %>%
  dplyr::mutate(
    declared_type_organism_name = na_if(declared_type_organism_name, ""),
    taxonomy_check_method = "NCBI_ANI"
  ) %>% 
  dplyr::arrange(desc(taxonomy_check_status), species_name) %>% 
  dplyr::select(-synonym_GB, -organism_name)

# all(ncbiTaxCheck$sampleId %in% genomeMetadata$sampleId)

## use inhouse data to perform taxonomy checking: because of additional type strains available
inhouseDf <- dplyr::filter(
  genomeMetadata, taxonomy_check_status != "OK" | 
    (source != "NCBI" & is.na(taxonomy_check_status))) %>% 
  dplyr::select(sampleId, source, AssemblyAccession, SpeciesName) %>% 
  dplyr::left_join(
    y = typeStrains, by = c("SpeciesName" = "type_organism_name")
  ) %>% 
  dplyr::left_join(
    y = typeStrainAni, by = c("sampleId" = "id1", "typeStrainId")
  ) %>% 
  # dplyr::filter(!is.na(ani) | is.na(typeStrainId)) %>%
  dplyr::group_by(sampleId) %>% 
  dplyr::arrange(dplyr::desc(ani), .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::select(
    sampleId, source, AssemblyAccession, SpeciesName, declared_type_id = typeStrainId, 
    declared_type_ANI = ani, declared_type_organism_name = type_organism_name
  )

if(!all(ncbiTaxCheck$sampleId %in% inhouseDf$sampleId)){
  stop(setdiff(ncbiTaxCheck$sampleId, inhouseDf$sampleId))
}

## get best ANI with its respective type species name
inhouseTaxCheck <- dplyr::left_join(
  inhouseDf, y = typeStrainAni, by = c("sampleId" = "id1")
) %>% 
  dplyr::group_by(sampleId) %>% 
  dplyr::arrange(dplyr::desc(ani), .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup() %>% 
  dplyr::rename(
    species_name = SpeciesName,
    best_match_species_name = type_organism_name,
    best_match_type_id = typeStrainId,
    best_match_type_ANI = ani
  ) %>% 
  dplyr::mutate(
    taxonomy_check_status = dplyr::case_when(
      declared_type_organism_name == best_match_species_name & 
        best_match_type_ANI >= 95 ~ "OK",
      declared_type_organism_name == best_match_species_name & 
        best_match_type_ANI < 95 ~ "Inconclusive",
      is.na(declared_type_id) ~ "Inconclusive",
      TRUE ~ "Failed"
    ),
    best_match_status = dplyr::case_when(
      declared_type_organism_name == best_match_species_name & 
        best_match_type_ANI < 96 ~ "below-threshold-match",
      declared_type_organism_name == best_match_species_name & 
        best_match_type_ANI >= 96 ~ "match",
      TRUE ~ "mismatch"
    ),
    comment = dplyr::case_when(
      is.na(declared_type_organism_name) ~ "no type strain to match",
      declared_type_ANI == 100 ~ "this is a type strain"
    ),
    taxonomy_check_method = "inhouse_fastANI",
    source = forcats::fct_relevel(source, "NCBI")
  ) %>% 
  dplyr::select(-mapped, -total, -best_match_type_id, -declared_type_id)



taxCheckDf <- dplyr::bind_rows(ncbiTaxCheck, inhouseTaxCheck) %>% 
  dplyr::group_by(sampleId) %>%
  dplyr::arrange(dplyr::desc(best_match_type_ANI), .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::rename(taxonomy_check_status_inhouse = taxonomy_check_status) %>% 
  dplyr::arrange(source, sampleId) %>% 
  dplyr::mutate(
    taxonomy_corrected = NA,
    new_species_name = NA
  ) %>% 
  dplyr::relocate(
    sampleId, taxonomy_check_status_inhouse,taxonomy_check_method,
    .before = taxonomy_corrected
  )

readr::write_tsv(x = taxCheckDf, file = confs$analysis$qc$files$tax_check)

################################################################################



