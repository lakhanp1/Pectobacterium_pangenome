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
  
  
  virSummary <- suppressMessages(
    readr::read_tsv(
      file = paste(confs$data$genomad$dir, "/", sampleId, "/", sampleId,
                   "_summary/", sampleId, "_virus_summary.tsv", sep = "")
    )
  ) %>% 
    dplyr::select(-n_genes) %>% 
    dplyr::rename_with(
      .fn = ~paste("genomad.", .x, sep = ""),
      .cols = c(virus_score, fdr, n_hallmarks, marker_enrichment)
    ) %>% 
    dplyr::mutate(
      sampleId = !!sampleId,
      SpeciesName = !!sampleInfoList[[sampleId]]$SpeciesName,
      chr = stringr::str_replace(seq_name, "(.*)\\|.*", "\\1"),
      .before = coordinates
    ) %>% 
    tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-")
  
  if(nrow(virSummary) > 0){
    checkv <- suppressMessages(
      readr::read_tsv(
        file = paste(confs$data$genomad$dir, "/", sampleId, "/", sampleId,
                     "_checkv/quality_summary.tsv", sep = "")
      )
    ) %>% 
      dplyr::left_join(y = virSummary, by = c("contig_id" = "seq_name"))
    
    prophageDf <- dplyr::bind_rows(prophageDf, checkv)
  }
  
}

dplyr::select(
  prophageDf, contig_id, sampleId, SpeciesName, chr, start, end, topology, taxonomy, 
  starts_with("genomad."), everything()
) %>% 
  # dplyr::glimpse() %>% 
  readr::write_tsv(file = confs$data$genomad$files$prophages)




