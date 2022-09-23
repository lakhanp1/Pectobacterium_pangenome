suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(XML))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(spData))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## Script to parse NCBI assembly and BioSample XML file
## summarize raw data and generate tables

rm(list = ls())

data(world, package = "spData")

#####################################################################
analysisName <- "raw_data_summary"
outDir <- here::here("analysis", "02_raw_data_summary")

file_assembly <- here::here("data/reference_data", "assembly_docsum.xml")
file_biosample <- here::here("data/reference_data", "biosample_docsum.xml")
file_inhouseMmetadata <- here::here("data/reference_data", "inhouse_samples_metadata.txt")
file_buscopMqc <- here::here("analysis/01_multiqc", "busco_prot_multiqc_data/multiqc_busco.txt")
file_buscogMqc <- here::here("analysis/01_multiqc", "busco_geno_multiqc_data/multiqc_busco.txt")
file_quastMqc <- here::here("analysis/01_multiqc", "quast_multiqc_data/multiqc_quast.txt")
file_ncbiAni <- here::here("data/other", "ANI_report_prokaryotes.txt")

#####################################################################
## inhouse genomes and QC data
inhouseGenomes <- suppressMessages(readr::read_tsv(file = file_inhouseMmetadata)) %>% 
  dplyr::mutate(
    source = "inhouse",
    SubmissionDate = "2023"
  )

buscopMqc <- suppressMessages(readr::read_tsv(file = file_buscopMqc)) %>% 
  dplyr::mutate(
    sampleId = stringr::str_replace(
      string = Sample, pattern = "short_summary.specific.enterobacterales_odb10.",
      replacement = ""
    )
  ) %>% 
  dplyr::select(-Sample, -lineage_dataset) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = c(starts_with("complete"), missing, fragmented),
      .fns = ~ round(. * 100/total, digits = 3)
    )
  ) %>% 
  dplyr::rename_with(
    .fn = ~ paste("buscop.", .x, sep = ""),
    .cols = !sampleId
  )

buscogMqc <- suppressMessages(readr::read_tsv(file = file_buscogMqc)) %>% 
  dplyr::mutate(
    sampleId = stringr::str_replace(
      string = Sample, pattern = "short_summary.specific.enterobacterales_odb10.",
      replacement = ""
    )
  ) %>% 
  dplyr::select(-Sample, -lineage_dataset) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = c(starts_with("complete"), missing, fragmented),
      .fns = ~ round(. * 100/total, digits = 3)
    )
  ) %>% 
  dplyr::rename_with(
    .fn = ~ paste("buscog.", .x, sep = ""),
    .cols = !sampleId
  )

quastCols <- c(
  length = "Total length",
  GC_per = "GC (%)",
  n_contigs = "# contigs",
  N50 = "N50",
  N90 = "N90",
  L50 = "L50",
  L90 = "L90",
  auN = "auN",
  N_per_100kb = "# N's per 100 kbp",
  contigs_gt_0kb = "# contigs (>= 0 bp)",
  contigs_gt_1kb = "# contigs (>= 1000 bp)",
  contigs_gt_5kb = "# contigs (>= 5000 bp)",
  contigs_gt_10kb = "# contigs (>= 10000 bp)",
  contigs_gt_25kb = "# contigs (>= 25000 bp)",
  contigs_gt_50kb = "# contigs (>= 50000 bp)",
  length_gt_0kb = "Total length (>= 0 bp)",
  length_gt_1kb = "Total length (>= 1000 bp)",
  length_gt_5kb = "Total length (>= 5000 bp)",
  length_gt_10kb = "Total length (>= 10000 bp)",
  length_gt_25kb = "Total length (>= 25000 bp)",
  length_gt_50kb = "Total length (>= 50000 bp)",
  largest_contig = "Largest contig"
)

quastMqc <- suppressMessages(readr::read_tsv(file = file_quastMqc)) %>% 
  dplyr::select(sampleId = Sample, !!!quastCols)

ncbiAni <- suppressMessages(readr::read_tsv(file = file_ncbiAni))

#####################################################################
## extract assembly information
asmDoc <- XML::xmlParse(file = file_assembly)
## asmDoc <- XML::xmlTreeParse(file = file_assembly)
asmNodes <- XML::getNodeSet(doc = asmDoc, path = "//DocumentSummary")
asmDf <- XML::xmlToDataFrame(nodes = asmNodes)

# xmlList <- XML::xmlToList(node = asmNodes[[1]])
# intersect(unlist(xmlList$PropertyList), c("latest_genbank", "latest_refseq"))

parse_assembly_docsum <- function(node){
  xmlList <- XML::xmlToList(node = node)
  
  exclFromRefSeq <- NA
  if(!is.null(node[["ExclFromRefSeq"]])){
    exclFromRefSeq <- paste(
      XML::xmlSApply(node[["ExclFromRefSeq"]], xmlValue),
      sep = ";", collapse = ";"
    )
    warning(xmlList$AssemblyAccession, exclFromRefSeq, ": ExclFromRefSeq\n")
  }
  
  Anomalous <- NA
  if(!is.null(node[["AnomalousList"]])){
    Anomalous <- paste(
      XML::xmlSApply(node[["AnomalousList"]][["Anomalous"]], xmlValue),
      sep = ";", collapse = ";"
    )
    warning(xmlList$AssemblyAccession, Anomalous, ": Anomalous\n")
  }
  
  replaced <- NA
  if(is.element(el = "replaced", xmlSApply(node[["PropertyList"]], xmlValue))){
    replaced <- "replaced"
    warning(xmlList$AssemblyAccession, replaced, ": replaced\n")
  }
  
  
  
  return(c(
    Id = xmlList$Id,
    taxonomy_check_status = xmlList$Meta$`taxonomy-check-status`,
    representative_status = xmlList$Meta$`representative-status`,
    ExclFromRefSeq = exclFromRefSeq,
    Anomalous = Anomalous,
    replaced = replaced,
    synonym_GB = xmlList$Synonym$Genbank,
    synonym_RS = xmlList$Synonym$RefSeq,
    synonym_similarity = xmlList$Synonym$Similarity
  ))
  
}

asmMetadata <- purrr::map_dfr(
  .x = asmNodes,
  .f = ~ parse_assembly_docsum(node = .)
) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = where(is.character),
      .fns = ~ stringr::str_replace_all(
        string = .x,
        pattern = regex(
          pattern = "^\\s*(not.applicable|missing|unknown|not.provided|na)\\s*$",
          ignore_case = TRUE
        ),
        replacement = NA_character_
      )
    )
    
  )

#####################################################################
## extract BioSample information
bsDoc <- XML::xmlParse(file = file_biosample)
bsNodes <- XML::getNodeSet(doc = bsDoc, path = "//DocumentSummary")

colInfo <- NULL

# aNode <- bsNodes[[1]]
# aNodeList <- XML::xmlToList(node = aNode)

bsMetadata <- purrr::map_dfr(
  .x = bsNodes,
  .f = function(aNode){
    aNodeDoc <- XML::xmlParse(XML::toString.XMLNode(aNode))
    aNodeList <- XML::xmlToList(node = aNode)
    
    ## parse all the metadata from Attribute nodes
    attrDf <- purrr::map_dfr(
      .x = XML::getNodeSet(doc = aNodeDoc, path = "//Attribute"),
      .f = function(nd){
        return(c(
          XML::xmlAttrs(node = nd),
          value = XML::xmlValue(x = nd)
        ))
      }
    )
    
    ## cleanup attribute names which will be used as column names
    attrDf <- dplyr::mutate(
      .data = attrDf,
      attribute_name = tolower(
        stringr::str_replace_all(
          string = attribute_name, pattern = "(\\W+)", replacement = "_"
        )
      ),
      attribute_name = stringr::str_replace(
        string = attribute_name, pattern = "_$", replacement = ""
      ),
      harmonized_name = dplyr::coalesce(harmonized_name, attribute_name)
    )
    
    ## cleanup the attribute values
    attrDf <- dplyr::mutate(
      .data = attrDf,
      value = stringr::str_replace_all(
        string = value,
        pattern = regex(
          pattern = "^\\s*(not.applicable|missing|unknown|not.provided|na)\\s*$",
          ignore_case = TRUE
        ),
        replacement = NA_character_
      )
    )
    
    ## add the attribute name and description to global DF
    colInfo <<- dplyr::bind_rows(
      colInfo,
      dplyr::filter(attrDf, !is.na(value)) %>% 
        dplyr::select(attribute_name, harmonized_name, display_name)
    )
    
    return(c(
      BioSampleAccn = aNodeList$Accession,
      tibble::deframe(x = dplyr::select(attrDf, harmonized_name, value))
    ))
  }
)

colStats <- dplyr::count(colInfo, harmonized_name, display_name, sort = TRUE)

bsMetadata <- dplyr::select(
  bsMetadata, BioSampleAccn, colStats$harmonized_name
) %>% 
  tibble::add_column(geo_loc_country = NA, .after = "geo_loc_name") %>% 
  dplyr::mutate(
    geo_loc_country = stringr::str_replace(
      string = geo_loc_name, pattern = ":.*", replacement = ""
    )
  ) %>% 
  dplyr::mutate(
    geo_loc_country = dplyr::case_when(
      geo_loc_country == "USA" ~ "United States",
      geo_loc_country == "South Korea" ~ "Republic of Korea",
      geo_loc_country == "Russia" ~ "Russian Federation",
      geo_loc_country == "USSR" ~ "Russian Federation",
      geo_loc_country == "Yugoslavia" ~ "Serbia",
      TRUE ~ geo_loc_country
    )
  )


if(!all(is.element(na.omit(unique(bsMetadata$geo_loc_country)), world$name_long))){
  stop("Unmatched country name")
}


#####################################################################

keyColumns <- c(
  "AssemblyAccession", "AssemblyName", "SpeciesName", "SpeciesTaxid", "Organism", 
  "AssemblyStatus", "AssemblyType", "BioSampleAccn", "SubmissionDate", "SubmitterOrganization",
  "FtpPath_RefSeq", "FtpPath_GenBank",  "Id"
)

## merge all data
ncbiData <- dplyr::select(.data = asmDf, !!!keyColumns)  %>% 
  dplyr::mutate(
    AssemblyName = stringr::str_replace_all(
      string = AssemblyName, pattern = "\\s+", replacement = "_"
    ),
    source = "NCBI"
  ) %>% 
  tidyr::unite(col = sampleId, AssemblyAccession, AssemblyName, sep = "_", remove = F) %>% 
  dplyr::left_join(y = asmMetadata, by = "Id") %>%
  dplyr::left_join(y = bsMetadata, by = "BioSampleAccn")


genomeMetadata <- dplyr::bind_rows(
  ncbiData, inhouseGenomes
) %>% 
  dplyr::left_join(y = buscopMqc, by = "sampleId") %>% 
  dplyr::left_join(y = buscogMqc, by = "sampleId") %>% 
  dplyr::left_join(y = quastMqc, by = "sampleId") %>% 
  dplyr::mutate(
    AssemblyStatus = dplyr::case_when(
      source == "inhouse" & n_contigs == 1 ~ "Complete Genome",
      source == "inhouse" & n_contigs > 1 ~ "Contig",
      TRUE ~ AssemblyStatus
    )
  )

genomeMetadata <- dplyr::select(
  genomeMetadata,
  sampleId, source, !!!keyColumns, !!!colnames(quastMqc),
  starts_with("buscog."), starts_with("buscop."),
  everything()
)

genomeMetadata <- dplyr::mutate(
  .data = genomeMetadata,  
  SubmissionDate = lubridate::ymd_hm(SubmissionDate),
  submission_y = lubridate::year(SubmissionDate),
  submission_m = lubridate::month(SubmissionDate)
)

#####################################################################
## taxonomy check failed/inconclusive information
taxCheckFail <- dplyr::select(
  ncbiData, sampleId, AssemblyAccession, synonym_GB, taxonomy_check_status
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
  dplyr::arrange(desc(taxonomy_check_status), species_name)

#####################################################################
readr::write_tsv(
  x = genomeMetadata,
  file = here::here("data/reference_data", "sample_metadata.tsv")
)

readr::write_tsv(
  x = taxCheckFail,
  file = file.path(outDir, "ncbi_taxanomy_qc_failed.tsv")
)

readr::write_tsv(
  x = colStats,
  file = file.path(outDir, "metadata_fields_stats.tsv")
)

wb <- openxlsx::createWorkbook()

currentSheet <- "metadata"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = genomeMetadata, withFilter = TRUE,
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

currentSheet <- "taxonomy_failed"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeData(
  wb = wb, sheet = currentSheet,
  x = stringr::str_c("#details for genomes where taxonomy check by NCBI failed")
)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = taxCheckFail, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA", startRow = 2
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 3, firstActiveCol = 2)


# openxlsx::openXL(wb)
saveWorkbook(
  wb = wb, file = here::here("data/reference_data", "sample_metadata.xlsx"), overwrite = TRUE
)

#####################################################################





