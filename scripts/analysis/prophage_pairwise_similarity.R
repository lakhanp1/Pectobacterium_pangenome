#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(configr))

# compute all-vs-all syntenic Jaccard index for consolidated prophages

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/compare_hg_sets.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
prophageLenCutoff <- confs$analysis$prophages$cutoff_length
maxSyntenyGap <- 2
minSyntenyChain <- 5

outDir <- confs$analysis$prophages$preprocessing$dir

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

logger::log_threshold(DEBUG)
file_log <- paste(
  "logs/prophage_similarity.", format(Sys.time(), "%Y%m%d.%H%M"), ".log",
  sep = ""
)
logger::log_appender(
  appender = logger::appender_file(file = file_log)
)

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$raw_prophage_hg)
) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName, genomeId), by = "sampleId"
  )

phages <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$consolidated_phages)
) %>% 
  dplyr::filter(filtered == 0) %>% 
  dplyr::mutate(
    fragments = stringr::str_split(fragments, ";")
  ) %>% 
  dplyr::select(
    prophage_id, fragments, nFragments, prophage_length, nHg, genomeId,
    parent, parentGenome
  )

phageList <- purrr::transpose(phages) %>% 
  purrr::set_names(nm = purrr::map(., "prophage_id"))

proHgL <- purrr::transpose(proHgs) %>%
  purrr::set_names(nm = purrr::map(., "prophage_id"))

################################################################################
# all vs all comparison

tmpFile <- paste(
  confs$analysis$prophages$preprocessing$files$phage_relationships,
  ".tmp", sep = ""
)

for (i in 1:nrow(phages)) {
  
  phageCmpDf <- tibble::tibble()
  
  for (j in i:nrow(phages)) {
    
    if (i == j) {
      next
    }
    
    phagePair <- c(phages$prophage_id[i], phages$prophage_id[j])
    
    logger::log_debug(
      "prophages[", i, ", ", j, "] = ", phagePair[1], ", ", phagePair[2]
    )
    
    # phagePair <- c("g_84.merged_vir_1", "g_68.vir_2")
    p1 <- phageList[[phagePair[1]]]
    p2 <- phageList[[phagePair[2]]]
    
    
    p1Hgs <- purrr::map(proHgL[p1$fragments], "hgs") %>% unlist(use.names = FALSE)
    p1nHgs <- length(p1Hgs)
    
    p2Hgs <- purrr::map(proHgL[p2$fragments], "hgs") %>% unlist(use.names = FALSE)
    p2nHgs <- length(p2Hgs)
    
    sharedHgs <- intersect(p1Hgs, p2Hgs)
    
    unionHgs <- union(p1Hgs, p2Hgs)
    
    # handle union and intersection calculation for prophages that have duplicated HGs
    if (any(duplicated(p1Hgs)) || any(duplicated(p2Hgs))) {
      # p1Hgs <- c("a", "b", "b", "c", "d", "d", "e", "f", "f", "g")
      # p2Hgs <- c("c", "d", "d", "e", "f", "f", "g", "g", "h", "i", "i", "j")
      p1Freq <- table(p1Hgs)
      p2Freq <- table(p2Hgs)
      
      # duplicated shared HGs
      sharedDups <- intersect(names(which(p1Freq > 1)), names(which(p2Freq > 1)))
      sharedDupFreq <- pmin(p1Freq[sharedDups], p2Freq[sharedDups]) - 1
      
      if (length(sharedDupFreq) > 0) {
        dupSharedHgs <- rep(names(sharedDupFreq), sharedDupFreq)
        
        sharedHgs <- append(sharedHgs, dupSharedHgs)
        unionHgs <- append(unionHgs, dupSharedHgs)
      }
      
      # duplicated unique HGs
      uniqueDupFreq <- c(p1Freq[setdiff(p1Hgs, sharedHgs)], p2Freq[setdiff(p2Hgs, sharedHgs)]) - 1
      
      if (length(uniqueDupFreq) > 0) {
        dupUniqueHgs <- rep(names(uniqueDupFreq), uniqueDupFreq)
        
        unionHgs <- append(unionHgs, dupUniqueHgs)
      }
      
      
    }
    
    
    # calculate syntenic HGs
    pairSyn <- list()
    
    if(length(sharedHgs) > 0){
      
      for (p1f in p1$fragments) {
        for (p2f in p2$fragments) {
          
          # cat(p1f, ":", proHgL[[p1f]]$hgs, "\n")
          # cat(p2f, ":", proHgL[[p2f]]$hgs, "\n\n")
          
          # get the shared syntenic homology groups between child and parent
          synteny <- syntenic_hg_overlap(
            ref = proHgL[[p1f]]$hgs,
            qur = proHgL[[p2f]]$hgs,
            minChainLen = min(proHgL[[p1f]]$nHgs, proHgL[[p2f]]$nHgs, minSyntenyChain),
            maxGapLen = maxSyntenyGap
          )
          
          if(!is.null(synteny)){
            pairSyn <- append(pairSyn, list(synteny))
          }
          
        }
      }
      
    }
    
    syntenicShared <- purrr::map(pairSyn, "lcs") %>%
      unlist()
    
    thisCmp <- list(
      phage1 = p1$prophage_id,
      phage2 = p2$prophage_id,
      phage1nHgs = p1nHgs,
      phage2nHgs = p2nHgs,
      nSharedHgs = length(sharedHgs),
      nSyntenicHgs = length(syntenicShared)
    )
    
    thisCmp$phage1Shared <- round(thisCmp$nSyntenicHgs / p1$nHg, digits = 4)
    thisCmp$phage2Shared <- round(thisCmp$nSyntenicHgs / p2$nHg, digits = 4)
    thisCmp$jaccardIndex <- round(
      thisCmp$nSyntenicHgs / length(unionHgs), digits = 4
    )
    thisCmp$contentDissimilarity <- round(
      1 - mean(c(thisCmp$phage1Shared, thisCmp$phage2Shared)), digits = 4
    )
    
    phageCmpDf <- dplyr::bind_rows(
      phageCmpDf,
      tibble::as_tibble(thisCmp)
    )
    
  }
  
  readr::write_tsv(
    x = phageCmpDf, file = tmpFile,
    col_names = (i == 1),
    append = (i > 1)
  )
  
}

file.copy(
  tmpFile, confs$analysis$prophages$preprocessing$files$phage_relationships,
  overwrite = TRUE
)






