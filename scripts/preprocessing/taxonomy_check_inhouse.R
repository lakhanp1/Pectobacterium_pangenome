suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(openxlsx))

## check the taxonomy for bacterial genomes
## generate a table for manual inspection

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- confs$analysis$qc$path

cutoff_ani_species <- 95L
cutoff_ani_ncbi <- 96L
################################################################################

genomeMetadata <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$prebuild_metadata)
) %>%
  dplyr::filter(
    ## no need to include assemblies that already have problem
    dplyr::if_all(
      .cols = c(ExclFromRefSeq, Anomalous, replaced),
      .fns = ~ is.na(.x)
    )
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
      .fns = ~ stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1")
    )
  )

typeStrainAni <- dplyr::left_join(
  typeStrains, aniDf,
  by = c("typeStrainId" = "id2")
)

################################################################################
## taxonomy check failed/inconclusive information

## NCBI data
ncbiTaxCheck <- dplyr::filter(genomeMetadata, source == "NCBI") %>%
  dplyr::select(
    sampleId, synonym_GB,
    ncbi_tax_status = taxonomy_check_status
  ) %>%
  dplyr::left_join(
    y = dplyr::select(
      ncbiAni, genbank_accession,
      ncbi_declared_species = declared_type_organism_name,
      ncbi_declared_ANI = declared_type_ANI,
      ncbi_declared_qcov = declared_type_qcoverage,
      ncbi_declared_scov = declared_type_scoverage,
      ncbi_best_ts_species = best_match_species_name,
      ncbi_best_ts_ANI = best_match_type_ANI,
      ncbi_best_ts_qcov = best_match_type_qcoverage,
      ncbi_best_ts_scov = best_match_type_scoverage,
      ncbi_best_ts_status = best_match_status,
      comment_ncbi = comment
    ),
    by = c("synonym_GB" = "genbank_accession")
  ) %>%
  dplyr::mutate(
    ncbi_declared_species = na_if(ncbi_declared_species, ""),
    ncbi_taxonomy_check_method = "NCBI_ANI"
  ) %>%
  dplyr::select(-synonym_GB)

################################################################################
## internal taxonomy checking
# all(ncbiTaxCheck$sampleId %in% genomeMetadata$sampleId)

## use inhouse data to perform taxonomy checking: because of additional type strains available
inhouseTypeAni <- dplyr::select(genomeMetadata, sampleId, AssemblyAccession, source, SpeciesName) %>%
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
    sampleId, source, AssemblyAccession, SpeciesName,
    inhouse_declared_id = typeStrainId,
    inhouse_declared_ANI = ani,
    inhouse_declared_species = type_organism_name
  )

if (!all(ncbiTaxCheck$sampleId %in% inhouseTypeAni$sampleId)) {
  stop(setdiff(ncbiTaxCheck$sampleId, inhouseTypeAni$sampleId))
}

## get best ANI with its respective type species name
inhouseTaxCheck <- dplyr::left_join(
  inhouseTypeAni,
  y = typeStrainAni, by = c("sampleId" = "id1")
) %>%
  dplyr::group_by(sampleId) %>%
  dplyr::arrange(dplyr::desc(ani), .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::rename(
    inhouse_best_ts_species = type_organism_name,
    inhouse_best_ts_id = typeStrainId,
    inhouse_best_ts_ANI = ani
  ) %>%
  dplyr::mutate(
    inhouse_tax_status = dplyr::case_when(
      inhouse_declared_species == inhouse_best_ts_species &
        inhouse_best_ts_ANI >= cutoff_ani_species ~ "OK",
      inhouse_declared_species == inhouse_best_ts_species &
        inhouse_best_ts_ANI < cutoff_ani_species ~ "Inconclusive",
      is.na(inhouse_declared_id) ~ "Inconclusive",
      TRUE ~ "Failed"
    ),
    inhouse_match_status = dplyr::case_when(
      inhouse_declared_species == inhouse_best_ts_species &
        inhouse_best_ts_ANI >= cutoff_ani_species ~ "match",
      inhouse_declared_species == inhouse_best_ts_species &
        inhouse_best_ts_ANI < cutoff_ani_species ~ "below-threshold-match",
      TRUE ~ "mismatch"
    ),
    comment_inhouse = dplyr::case_when(
      is.na(inhouse_declared_species) ~ "no type strain to match",
      inhouse_declared_ANI == 100 ~ "this is a type strain"
    ),
    inhouse_taxonomy_check_method = "inhouse_fastANI"
  ) %>%
  dplyr::select(-mapped, -total, -inhouse_best_ts_id, -inhouse_declared_id)

taxCheckDf <- dplyr::left_join(
  inhouseTaxCheck, ncbiTaxCheck,
  by = "sampleId"
) %>%
  dplyr::filter(
    inhouse_tax_status != "OK" | ncbi_tax_status != "OK" | source != "NCBI"
  ) %>%
  tidyr::unite(col = "comment", comment_inhouse, comment_ncbi) %>%
  dplyr::mutate(
    ncbi_eq_inhouse = if_else(
      ncbi_best_ts_species == inhouse_best_ts_species, true = "TRUE", false = "FALSE"
    )
  )

bestAni <- dplyr::select(taxCheckDf, c(sampleId, ends_with("_ANI"))) %>%
  tidyr::pivot_longer(
    cols = -sampleId,
    names_to = "best_ANI_from",
    values_to = "best_ANI",
    values_drop_na = TRUE
  ) %>%
  dplyr::group_by(sampleId) %>%
  dplyr::arrange(dplyr::desc(best_ANI), .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

taxCheckDf <- dplyr::left_join(
  taxCheckDf, bestAni,
  by = "sampleId"
) %>%
  dplyr::arrange(source, SpeciesName) %>%
  dplyr::mutate(
    best_ani_species = dplyr::case_when(
      best_ANI_from == "inhouse_declared_ANI" ~ inhouse_declared_species,
      best_ANI_from == "inhouse_best_ts_ANI" ~ inhouse_best_ts_species,
      best_ANI_from == "ncbi_declared_ANI" ~ ncbi_declared_species,
      best_ANI_from == "ncbi_best_ts_ANI" ~ ncbi_best_ts_species,
      TRUE ~ "NA"
    ),
    ani_matches = dplyr::if_else(
      SpeciesName == best_ani_species, true = TRUE, false = FALSE
    ),
    taxonomy_check_method = dplyr::case_when(
      best_ANI >= cutoff_ani_species &
        (best_ANI_from %in% c("inhouse_declared_ANI", "inhouse_best_ts_ANI")) ~
        inhouse_taxonomy_check_method,
      best_ANI >= cutoff_ani_species &
        (best_ANI_from %in% c("ncbi_declared_ANI", "ncbi_best_ts_ANI")) ~
        ncbi_taxonomy_check_method,
      TRUE ~ ""
    ),
    taxonomy_check_status_inhouse = dplyr::case_when(
      best_ANI_from %in% c("inhouse_declared_ANI", "inhouse_best_ts_ANI") ~
        inhouse_tax_status,
      best_ANI_from %in% c("ncbi_declared_ANI", "ncbi_best_ts_ANI") ~
        ncbi_tax_status
    )
  ) %>%
  dplyr::mutate(
    taxonomy_corrected = NA,
    new_species_name = NA
  )

glimpse(taxCheckDf)


################################################################################
readr::write_tsv(x = taxCheckDf, file = confs$analysis$qc$files$tax_check)

file_xls <- paste(
  tools::file_path_sans_ext(confs$analysis$qc$files$tax_check), ".xlsx",
  sep = ""
)

negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
neutralStyle <- openxlsx::createStyle(fontColour = "#0000ff", bgFill = "#ccccff")
naStyle <- openxlsx::createStyle(fontColour = "black", bgFill = "#e6e6e6")

currentSheet <- "tax_check"

if (file.exists(file_xls)) {
  wb <- openxlsx::loadWorkbook(file = file_xls)
  removeWorksheet(wb, sheet = currentSheet)
} else {
  wb <- openxlsx::createWorkbook()
}
openxlsx::addWorksheet(wb, sheetName = currentSheet)

## tax_status column formating

openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = taxCheckDf, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)

for (col in which(colnames(taxCheckDf) %in% c("inhouse_tax_status", "ncbi_tax_status"))) {
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "contains", rule = "OK", style = posStyle
  )
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "contains", rule = "Failed", style = negStyle
  )
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "contains", rule = "Inconclusive", style = neutralStyle
  )
}

for (col in which(
  colnames(taxCheckDf) %in%
    c(
      "inhouse_declared_ANI", "inhouse_best_ts_ANI", "ncbi_declared_ANI",
      "ncbi_best_ts_ANI", "best_ANI"
    )
)) {
  ## ANI column formatting
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "expression", rule = ">=95", style = posStyle
  )
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "expression", rule = "<95", style = negStyle
  )
  openxlsx::conditionalFormatting(
    wb = wb, sheet = currentSheet,
    cols = col, rows = 1:nrow(taxCheckDf) + 1,
    type = "contains", rule = "NA", style = naStyle
  )
}


openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 4)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = file_xls, overwrite = TRUE
)

################################################################################
