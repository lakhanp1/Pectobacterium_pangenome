suppressMessages(library(tidyverse))
suppressMessages(library(Rsamtools))

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

# sampleId <- "JGAR-100719-1"
prophageDf <- NULL

for (sampleId in sampleInfo$sampleId) {
  
  virPathPrefix <- paste(confs$data$genomad$dir, "/", sampleId, "/", sampleId, sep = "")
  
  virSummary <- suppressMessages(
    readr::read_tsv(
      file = paste(virPathPrefix, "_summary/", sampleId, "_virus_summary.tsv", sep = "")
    )
  ) %>% 
    dplyr::select(-n_genes) %>% 
    dplyr::rename_with(
      .fn = ~paste("genomad.", .x, sep = ""),
      .cols = c(virus_score, fdr, n_hallmarks, marker_enrichment)
    ) %>% 
    dplyr::mutate(
      sampleId = !!sampleId,
      Genome = !!sampleInfoList[[sampleId]]$Genome,
      SpeciesName = !!sampleInfoList[[sampleId]]$SpeciesName,
      chr = stringr::str_replace(seq_name, "(.*)\\|.*", "\\1"),
      .before = coordinates
    ) %>% 
    tidyr::separate(col = coordinates, into = c("start", "end"), sep = "-")
  
  if(nrow(virSummary) > 0){
    checkv <- suppressMessages(
      readr::read_tsv(
        file = paste(virPathPrefix, "_checkv/quality_summary.tsv", sep = "")
      )
    ) %>% 
      dplyr::left_join(y = virSummary, by = c("contig_id" = "seq_name")) %>% 
      dplyr::mutate(
        prophage_id = paste("g_", Genome, ".vir_", 1:n(), sep = "")
      )
    
    prophageDf <- dplyr::bind_rows(prophageDf, checkv)
    
  }
  
}

prophageDf <- dplyr::select(
  prophageDf, contig_id, prophage_id, sampleId, SpeciesName, chr, start, end,
  topology, taxonomy, starts_with("genomad."), everything()
)

readr::write_tsv(prophageDf, file = confs$data$genomad$files$prophages)

dplyr::select(prophageDf, sampleId, chr, start, end, prophage_id) %>% 
  dplyr::mutate(
    fna = paste(confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""),
    region = paste(chr, ":", start, "-", end, sep = ""),
    out = paste(confs$data$genomad$dir, "/phage_seqs/", prophage_id, ".fna", sep = ""),
    faidx = paste("samtools faidx", fna, region, ">", out)
  ) %>% 
  dplyr::select(faidx) %>% 
  readr::write_tsv(
    file = "./scripts/a_preprocessing/prophage_genomes_extract.sh",
    col_names = F
  )


