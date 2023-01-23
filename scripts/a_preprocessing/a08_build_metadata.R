suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(openxlsx))

## prepare final metadata file by correcting taxonomy

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

cutoff_buscog <- confs$parameters$cutoff_busco

################################################################################

genomeMetadata <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$prebuild_metadata)
)

correctTaxo <- suppressMessages(
  readr::read_tsv(confs$data$reference_data$files$taxo_correction)
) %>% 
  dplyr::mutate(
    taxonomy_check_status_inhouse = dplyr::if_else(
      condition = taxonomy_corrected == TRUE,
      true = "Corrected", false = taxonomy_check_status_inhouse
    )
  )

duplicateGenomes <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$duplicate_genomes)
) %>%
  dplyr::filter(identical == TRUE)
excludeGenomes <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$exclude_genomes, col_names = "sampleId")
)

################################################################################
correctMeta <- dplyr::left_join(
  x = genomeMetadata, y = correctTaxo, by = "sampleId"
) %>% 
  dplyr::mutate(
    taxonomy_check_status = dplyr::if_else(
      condition = !is.na(taxonomy_check_status_inhouse),
      true = taxonomy_check_status_inhouse, false = taxonomy_check_status
    ),
    SpeciesName = dplyr::if_else(
      condition = !is.na(new_species_name),
      true = new_species_name, false = SpeciesName
    )
  ) %>% 
  tidyr::replace_na(
    replace = list(taxonomy_corrected = FALSE, taxonomy_check_method = "NCBI_ANI")
  ) %>% 
  dplyr::select(-taxonomy_check_status_inhouse, -new_species_name) %>% 
  dplyr::relocate(taxonomy_check_method, .after = taxonomy_check_status) %>% 
  dplyr::relocate(taxonomy_corrected, .after = SpeciesName)



correctMeta %<>% 
  dplyr::mutate(
    filtered = dplyr::case_when(
      !is.na(ExclFromRefSeq) ~ "ExclFromRefSeq",
      !is.na(Anomalous) ~ "Anomalous",
      !is.na(replaced) ~ "replaced",
      sampleId %in% c(duplicateGenomes$genome2, excludeGenomes$sampleId) ~ "duplicate",
      buscog.complete < cutoff_buscog ~ paste("BUSCO.geno <", cutoff_buscog),
      TRUE ~ "PASS"
    ),
    type_material = stringr::str_replace(
      string = tolower(type_material), pattern = "(type (strain|material)).*", 
      replacement = "type strain"
    )
  ) %>% 
  dplyr::relocate(filtered, .after = AssemblyName)

#####################################################################
readr::write_tsv(
  x = correctMeta,
  file = confs$data$reference_data$files$metadata
)

## get column information from prebuild metadata
pbwb <- openxlsx::loadWorkbook(file = confs$analysis$qc$files$prebuild_metadata_xls)
colStats <- openxlsx::readWorkbook(pbwb, sheet = "column_info")

wb <- openxlsx::createWorkbook()

currentSheet <- "metadata"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = correctMeta, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

currentSheet <- "column_info"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = colStats, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = confs$data$reference_data$files$metadata_xls, overwrite = TRUE
)

################################################################################



