#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))

# correct predicted prophages using pangenome homology groups
# this will involve merging prophages in genome that were split because of 
# the assembly breaks

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/association_analysis.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)


treeMethod <- "kmer_nj"     #ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- confs$analysis$prophages$dir
outPrefix <- ""

orgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

################################################################################

prophageDf <- suppressMessages(readr::read_tsv(confs$data$genomad$files$prophages)) %>% 
  dplyr::select(prophage_id, prophage_length = length)

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>% 
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName, Genome, N50), by = "sampleId"
  ) %>% 
  dplyr::arrange(desc(nHgs), desc(N50)) %>%
  dplyr::left_join(y = prophageDf, by = "prophage_id") %>% 
  dplyr::filter(prophage_length >= 5000)

proHgL <-  purrr::transpose(proHgs) %>% 
  purrr::set_names(nm = purrr::map(., "prophage_id"))


sharedStats <- tibble::tibble() 

# ph <- proHgs$prophage_id[1]
# ph <- "g_99.vir_1"
# ph <- "g_261.vir_5"

for (ph in proHgs$prophage_id) {
  cat("finding child prophages for", ph, "\n")
  
  parentPh <- purrr::pluck(proHgL, !!!ph)
  
  foundChild <- FALSE
  proHgL[[ph]]$parent_tested <- TRUE
  
  # for each genome
  for (gn in sampleInfo$Genome) {
    
    if(parentPh$Genome == gn){
      next
    }
    
    gnPhs <- dplyr::filter(proHgs, Genome == !!gn) %>% 
      dplyr::pull(prophage_id)
    
    childPhages <- NULL
    
    # for each prophage in genome, check for overlap
    for (gp in gnPhs) {
      
      # skip if this prophage was used as parent before
      if(isTRUE(proHgL[[gp]]$parent_tested)){
        next
      }
      
      gpL <- purrr::pluck(proHgL, !!!gp)
      
      nUniqHgs <- length(setdiff(gpL$hgs, parentPh$hgs))
      nSharedHgs <- length(intersect(gpL$hgs, parentPh$hgs))
      
      # set parent
      if(nUniqHgs <= 2 & nSharedHgs > nUniqHgs){
        childPhages <- append(childPhages, gpL$prophage_id)
      }
      
    }
    
    # further checks and filtering on child phages
    if(!is.null(childPhages)){
      
      childHgs <- purrr::map(proHgL[childPhages], "hgs") %>% unlist() %>% unique()
      nSharedHgs <- length(intersect(childHgs, parentPh$hgs))
      nUniqParent <- length(setdiff(parentPh$hgs, childHgs))
      nUniqChild <- length(setdiff(childHgs, parentPh$hgs))
      childPhageLen <- purrr::map_dbl(proHgL[childPhages], "prophage_length") %>% sum()
      
      contigCheck <- TRUE
      
      # check that all phages are located on independent contigs or at 5'/3' 
      # ends of large contigs
      if(length(childPhages) > 1){
        
        childLoc <- get_hg_sets_location(
          hgs = proHgL[["g_378.vir_35"]]$hgs, genome = gn, pandb = orgDb
        )
        
        phageContigPos <- purrr::map_dbl(
          proHgL[childPhages],
          .f = ~ get_hg_sets_location(hgs = .x$hgs, genome = gn, pandb = orgDb)
        ) 
          
        if(any(phageContigPos == -1)){
          contigCheck <- FALSE
        }
      }
      
      relation <- dplyr::case_when(
        nSharedHgs / parentPh$nHgs == 1 ~ "identical",
        nSharedHgs / parentPh$nHgs >= 0.9 ~ "high_quality",
        nSharedHgs / parentPh$nHgs >= 0.5 ~ "medium_quality",
        nSharedHgs / parentPh$nHgs < 0.5 ~ "low_quality",
        contigCheck == FALSE ~ "NA",
        TRUE ~ "NA"
      )
      
      if(relation != "NA"){
        foundChild <- TRUE
        
        gnPhageRel <- tibble::tibble(
          parent = parentPh$prophage_id,
          parentGenome = parentPh$Genome,
          nHgParent = parentPh$nHgs,
          hasChild = foundChild,
          relation = relation,
          child = paste(childPhages, collapse = ";"),
          childGenome = gn,
          nChildPhages = length(childPhages),
          nHgChild = length(childHgs),
          nShared = nSharedHgs,
          fracSharedParent = nShared/nHgParent,
          fracSharedChild = nShared/nHgChild
        )
        
        # print(gnPhageRel)
        
        # proHgL[[gp]]$parent <- append(proHgL[[gp]]$parent, parentPh$prophage_id)
        # proHgL[[ph]]$child <- append(proHgL[[ph]]$child, childPhages)
        
        sharedStats <- dplyr::bind_rows(sharedStats, gnPhageRel)
      }
      
    }
    
  }
  
  if(!foundChild){
    sharedStats <- dplyr::bind_rows(
      sharedStats,
      tibble::tibble(
        parent = parentPh$prophage_id,
        parentGenome = parentPh$Genome,
        nHgParent = parentPh$nHgs,
        hasChild = foundChild
      )
    )
  }
  
  
}


readr::write_tsv(
  x = sharedStats,
  file = file.path(outDir, "pangenome_phage_similarity.tab")
)



