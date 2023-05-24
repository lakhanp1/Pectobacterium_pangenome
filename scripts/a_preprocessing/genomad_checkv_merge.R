suppressMessages(library(tidyverse))

# process and combine genomad and checkv output files for all samples

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

sampleId <- "GCF_024498715.1_ASM2449871v1"
prophageDf <- NULL

for (sampleId in sampleInfo$sampleId) {
  checkv <- suppressMessages(
    readr::read_tsv(
      file = paste(confs$data$genomad$dir, "/", sampleId, "/", sampleId,
                   "_checkv/quality_summary.tsv", sep = "")
    )
  )
  
  virSummary <- suppressMessages(
    readr::read_tsv(
      file = paste(confs$data$genomad$dir, "/", sampleId, "/", sampleId,
                   "_summary/", sampleId, "_virus_summary.tsv", sep = "")
    )
  ) %>% 
    dplyr::select(seq_name, topology, coordinates) %>% 
    dplyr::mutate(
      sampleId = !!sampleId,
      chr = stringr::str_replace(seq_name, "(.*)\\|.*", "\\1"),
      .before = coordinates
    ) %>% 
    tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-") %>% 
    dplyr::left_join(y = checkv, by = c("seq_name" = "contig_id"))

  prophageDf <- dplyr::bind_rows(prophageDf, virSummary)
  
}

readr::write_tsv(x = prophageDf, file = confs$data$genomad$files$prophages)




