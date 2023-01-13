suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(taxize))
suppressPackageStartupMessages(library(XML))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(spData))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## Script to parse NCBI assembly and BioSample XML file
## summarize raw data and generate tables

rm(list = ls())

data(world, package = "spData")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
################################################################################
analysisName <- "raw_data_summary"

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- confs$analysis$qc$dir

################################################################################
## inhouse genomes and QC data
inhouseGenomes <- suppressMessages(
  readr::read_tsv(file = confs$data$reference_data$files$inhouse_metadata)
) %>% 
  dplyr::mutate(
    SubmissionDate = "2023",
    collection_date = as.character(collection_date)
  )

missingMetadata <- suppressMessages(
  readr::read_tsv(file = confs$data$other$files$pecto1_metadata, na = "-")
)

buscopMqc <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$busco_protein)
) %>% 
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

buscogMqc <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$busco_genome)
) %>% 
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

quastMqc <- suppressMessages(
  readr::read_tsv(file = confs$analysis$qc$files$quast)
) %>% 
  dplyr::select(sampleId = Sample, !!!quastCols)

ncbiAni <- suppressMessages(
  readr::read_tsv(file = confs$data$other$files$ani_report)
  )

################################################################################
## extract assembly information
asmDoc <- XML::xmlParse(file = confs$data$reference_data$files$assembly_xml)
## asmDoc <- XML::xmlTreeParse(file = confs$data$reference_data$files$assembly_xml)
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


################################################################################
## extract BioSample information
bsDoc <- XML::xmlParse(file = confs$data$reference_data$files$biosample_xml)
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
)

################################################################################

assemblyCols <- c(
  "AssemblyAccession", "AssemblyName", "SpeciesName", "SpeciesTaxid", "Organism", 
  "AssemblyStatus", "AssemblyType", "BioSampleAccn", "SubmissionDate", "SubmitterOrganization",
  "FtpPath_RefSeq", "FtpPath_GenBank",  "Id"
)

## merge all data
ncbiData <- dplyr::select(.data = asmDf, !!!assemblyCols)  %>% 
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
      source != "NCBI" & n_contigs == 1 ~ "Complete Genome",
      source != "NCBI" & n_contigs > 1 ~ "Contig",
      AssemblyStatus == "Chromosome" ~ "Complete Genome",
      TRUE ~ AssemblyStatus
    ),
    sampleName = dplyr::if_else(
      condition = is.na(AssemblyAccession), true = sampleId, false = AssemblyAccession
    )
  )


informativeCols <- c(
  "AssemblyAccession", "AssemblyName", "SpeciesName", "SpeciesTaxid", "Organism",
  "AssemblyStatus", "AssemblyType", "BioSampleAccn", "strain", "geo_loc_name",
  "collection_date", "host", "isolation_source",
  "env_broad_scale", "env_local_scale", "env_medium"
)

genomeMetadata <- dplyr::select(
  genomeMetadata,
  sampleId, sampleName, source, !!!informativeCols, !!!colnames(quastMqc),
  starts_with("buscog."), starts_with("buscop."),
  everything()
)


################################################################################
## add missing metadata for NCBI genomes: this metadata is collected
## from the people who submitted the genomes but did not include the metadata
missingMetadata <- dplyr::left_join(
  x = missingMetadata,
  y = dplyr::filter(genomeMetadata, source == "NCBI") %>% 
    dplyr::select(AssemblyAccession) %>% 
    dplyr::mutate(rowIndex = 1:n()),
  by = "AssemblyAccession"
)

for(mtcol in setdiff(colnames(missingMetadata), c("AssemblyAccession", "rowIndex"))){
  genomeMetadata[[mtcol]][missingMetadata$rowIndex] <- dplyr::coalesce(
    genomeMetadata[[mtcol]][missingMetadata$rowIndex],
    missingMetadata[[mtcol]]
  )
}


#####################################################################
## correct dates
genomeMetadata <- dplyr::mutate(
  .data = genomeMetadata,  
  SubmissionDate = lubridate::ymd_hm(SubmissionDate),
  submission_y = lubridate::year(SubmissionDate),
  collection_date = lubridate::ymd(collection_date, truncated = 2),
  collection_year = lubridate::year(collection_date)
) %>% 
  dplyr::relocate(collection_year, .before = collection_date) %>% 
  dplyr::relocate(submission_y, .before = SubmissionDate)

## correct country name
genomeMetadata %<>% tibble::add_column(geo_loc_country = NA, .after = "geo_loc_name") %>% 
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


if(!all(is.element(na.omit(unique(genomeMetadata$geo_loc_country)), world$name_long))){
  stop("Unmatched country name")
}


#####################################################################
readr::write_tsv(
  x = genomeMetadata,
  file = confs$analysis$qc$files$prebuild_metadata
)

## if exists the open the excel, else create new
wb <- openxlsx::createWorkbook()
if(file.exists(confs$analysis$qc$files$prebuild_metadata_xls)){
  wb <- openxlsx::loadWorkbook(file = confs$analysis$qc$files$prebuild_metadata_xls)
}

currentSheet <- "metadata"
if(currentSheet %in% names(wb)){
  openxlsx::removeWorksheet(wb, currentSheet)
}
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = genomeMetadata, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

currentSheet <- "column_info"
if(currentSheet %in% names(wb)){
  openxlsx::removeWorksheet(wb, currentSheet)
}
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = colStats, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
saveWorkbook(
  wb = wb,
  file = confs$analysis$qc$files$prebuild_metadata_xls, overwrite = TRUE
)

################################################################################






