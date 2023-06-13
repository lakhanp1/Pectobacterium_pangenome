#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(logger))

# build prophage hierarchy tree using the homology groups

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
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


logger::log_threshold(WARN)
logger::log_appender(
  appender = logger::appender_file(
    file = paste(outDir, "/prophage_correction.log", sep = "")
  )
)

################################################################################
# sorting by N50 to avoid the circular links. In case of duplicate prophages
# across multiple genome assemblies, prophage from assembly with highest N50 is
# always assigned as parent to the prophages from remaining assembly. Prophage
# from the highest quality assembly is root in this case
sampleInfo <- get_metadata(file = panConf$files$metadata) %>% 
  dplyr::arrange(N50)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

################################################################################

prophageDf <- suppressMessages(readr::read_tsv(confs$data$genomad$files$prophages)) %>% 
  dplyr::select(prophage_id, prophage_length = length, chr)

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
  dplyr::left_join(y = prophageDf, by = "prophage_id") %>% 
  dplyr::filter(prophage_length >= 5000) %>%
  dplyr::arrange(desc(nHgs), desc(N50))

proHgL <-  purrr::transpose(proHgs) %>% 
  purrr::set_names(nm = purrr::map(., "prophage_id"))


sharedStats <- tibble::tibble() 
childParentList <- list()

# sampleInfo <- dplyr::filter(
#   sampleInfo,
#   # Genome %in% "378",
#   Genome %in% c("99", "108", "109", "110", "111","112", "113", "114", "115")
# )
# gn <- "378"
# ph <- "g_111.vir_3"
# ph <- "g_261.vir_5"

for (gn in sampleInfo$Genome) {
  logger::log_info('processing phages in genome: ', gn)
  
  gnPhages <- dplyr::filter(proHgs, Genome == !!gn) %>% 
    dplyr::pull(prophage_id)
  
  # purrr::keep(.x = proHgL, .p = function(x) x$Genome == gn)
  
  bestParent <- list()
  
  # iterate over all phages to find the best parent phage
  for (ph in proHgs$prophage_id) {
    parentPh <- purrr::pluck(proHgL, !!!ph)
    
    if(parentPh$Genome == gn){
      next
    }
    
    childPhages <- NULL
    
    # for each prophage in genome, check for overlap with prophage list
    for (gp in gnPhages) {
      gpL <- purrr::pluck(proHgL, !!!gp)
      
      nUniqHgs <- length(setdiff(gpL$hgs, parentPh$hgs))
      nSharedHgs <- length(intersect(gpL$hgs, parentPh$hgs))
      
      if(nUniqHgs <= 2 && nSharedHgs > nUniqHgs){
        childPhages <- append(childPhages, gpL$prophage_id)
      }
      
    }
    
    # evaluate if this parent is the best, if not, update the parent
    if(!is.null(childPhages)){
      
      # prevent circular links: was current parentPh detected as child
      # to the current childPhages
      # need a recursion here
      if(!is.null(childParentList[[parentPh$prophage_id]])){
        if (setequal(childParentList[[parentPh$prophage_id]], childPhages)) {
          logger::log_warn('Circular link ', parentPh$prophage_id, ' and ',
                           paste(childPhages, collapse = "; "))
          next
        } 
      }
      
      logger::log_debug('parent phage: ', ph, ' | evaulating children: ',
                        paste(childPhages, collapse = "; "))
      
      childHgs = purrr::map(proHgL[childPhages], "hgs") %>% unlist() %>% unique()
      contigCheck <- TRUE
      
      # in case of multiple child-phages mapping to one parent prophage, check that
      # all phages are located on independent contigs or at 5'/3' ends of large contigs
      if(length(childPhages) > 1){
        
        phageContigPos <- purrr::map_dbl(
          proHgL[childPhages],
          .f = ~ get_hg_sets_location(hgs = .x$hgs, genome = gn, chr = .x$chr, pandb = orgDb)
        ) 
        
        if(any((phageContigPos != 0) & (phageContigPos != 1))){
          contigCheck <- FALSE
        }
        
      }
      
      # check for duplicate genome child-parent & parent-child circular links
      
      
      thisParent <- list(
        parent = parentPh$prophage_id,
        parentGenome = parentPh$Genome,
        nHgParent = parentPh$nHgs,
        children = childPhages,
        childGenome = gn,
        childString = paste(sort(childPhages), collapse = ";"),
        nChild = length(childPhages),
        childHgs = childHgs,
        nSharedHgs = length(intersect(childHgs, parentPh$hgs)),
        nUniqParent = length(setdiff(parentPh$hgs, childHgs)),
        nUniqChild = length(setdiff(childHgs, parentPh$hgs)),
        childPhageLen = purrr::map_dbl(proHgL[childPhages], "prophage_length") %>% sum()
      )
      
      thisParent$nHgChild <- length(thisParent$childHgs)
      thisParent$perSharedParent <- thisParent$nSharedHgs / thisParent$nHgParent
      thisParent$perSharedChild <- thisParent$nSharedHgs / thisParent$nHgChild
      
      thisParent$relation <- dplyr::case_when(
        thisParent$perSharedParent == 1 ~ "identical",
        thisParent$perSharedParent >= 0.9 ~ "high_quality",
        thisParent$perSharedParent >= 0.5 ~ "medium_quality",
        thisParent$perSharedParent < 0.5 ~ "low_quality",
        contigCheck == FALSE ~ "NA",
        TRUE ~ "NA"
      )
      
      logger::log_debug("contigCheck: ", contigCheck)
      
      # compare with previous best parent 
      if (contigCheck) {
        
        ## if no parent, make current phage as parent
        if (length(bestParent) == 0) {
          bestParent[1] <- list(thisParent)
          logger::log_info("First parent: ", thisParent$parent, " for children: ",
                           paste(thisParent$children, collapse = "; "))
          logger::log_info("#parents: ", length(bestParent))
          
          next
        }
        
        parentChildSetFound <- FALSE
        betterParentFound <- FALSE
        
        # iterate over parent set and replace the parent element combination 
        ## if the current one is better
        for (i in 1:length(bestParent)) {
          # make sure correct set is being compared
          if(any(thisParent$children %in% bestParent[[i]]$children) && !betterParentFound){
            
            parentChildSetFound <- TRUE
            
            oldParentPhage <- purrr::pluck(proHgL, !!!bestParent[[i]]$parent)
            
            compSharedHgs <- dplyr::case_when(
              thisParent$nSharedHgs == bestParent[[i]]$nSharedHgs ~ "equal",
              thisParent$nSharedHgs > bestParent[[i]]$nSharedHgs ~ "more",
              thisParent$nSharedHgs < bestParent[[i]]$nSharedHgs ~ "less"
            )
            
            compLen <- dplyr::case_when(
              thisParent$childPhageLen == bestParent[[i]]$childPhageLen ~ "equal",
              thisParent$childPhageLen > bestParent[[i]]$childPhageLen ~ "longer",
              thisParent$childPhageLen < bestParent[[i]]$childPhageLen ~ "smaller"
            )
            
            isBetterParent <- dplyr::case_when(
              compSharedHgs == "more" ~ "more_hgs",
              compSharedHgs == "equal" && compLen == "longer" ~ "longer_prophage",
              compSharedHgs == "equal" && parentPh$N50 > oldParentPhage$N50 ~ "better_N50",
              # compSharedHgs == "less" ~ FALSE,
              TRUE ~ "no"
            )
            
            if (isBetterParent != "no") {
              logger::log_info("Found better parent: sharedHgs = ", compSharedHgs,
                               ", length = ", compLen, " :", isBetterParent)
              logger::log_info(
                "nHgs ", thisParent$nSharedHgs, "==", bestParent[[i]]$nSharedHgs,
                ", Len ", thisParent$childPhageLen, "==", bestParent[[i]]$childPhageLen,
                ", N50: ", parentPh$N50, " > ", oldParentPhage$N50
              )
              logger::log_info(
                "Previous parent: #", i, ": ", bestParent[[i]]$parent, " for children: ",
                paste(bestParent[[i]]$children, collapse = ";")
              )
              logger::log_info(
                "New parent: #", i, ": ", thisParent$parent, " for children: ",
                paste(thisParent$children, collapse = ";")
              )
              
              betterParentFound <- TRUE
              bestParent[[i]] <- thisParent
              
              logger::log_info("#parents: ", length(bestParent))
              
            }
            
            # IMP to break the loop here
            # break
            
          } else if (betterParentFound) {
            # if better parent was found in one of the previous iteration,
            # ensure the other parent-child relationship that involve any
            # of the child in this iteration is set to NA 
            
            if(any(thisParent$children %in% bestParent[[i]]$children)){
              # cannot set to NULL as this will change the for loop index 
              logger::log_info(
                "Removing parent #", i, ": ", bestParent[[i]]$parent, " for children: ",
                paste(bestParent[[i]]$children, collapse = ";")
              )
              
              bestParent[[i]] <- NULL
              
              logger::log_info("#parents: ", length(bestParent))
              
            }
          }
          
        }
        
        # parent exists but not for this child phages, add the current phage to bestParentList
        if(!parentChildSetFound){
          bestParent <- append(bestParent, list(thisParent))
          
          logger::log_info("New parent: #", length(bestParent), ": ", thisParent$parent,
                           " for children: ", paste(thisParent$children, collapse = ";"))
          logger::log_info("#parents: ", length(bestParent))
        }
        
      }
      
    }
    
  }
  
  if (length(bestParent) > 0) {
    
    thisGnParents <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c("childString", "childGenome", "nHgChild", "parent", "parentGenome", "nHgParent",
        "nSharedHgs", "perSharedParent", "perSharedChild", "relation")
    ) %>% 
      dplyr::rename(child = childString)
    
    
    # for remaining phages where parent was not found
    noParents <- setdiff(gnPhages, unlist(purrr::map(bestParent, "children")))
    if (length(noParents)) {
      thisGnParents <- purrr::map_dfr(
        .x = proHgL[noParents], .f = `[`,
        c("prophage_id", "Genome", "nHgs")
      ) %>% 
        dplyr::rename(child = prophage_id, childGenome = Genome, nHgChild = nHgs) %>% 
        dplyr::bind_rows(thisGnParents)
    }
    
    # list for circular child-parent-child linking handling
    thisGnChildren <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c("children", "parent")
    ) %>% 
      tibble::deframe() %>% 
      as.list()
    
    childParentList <- append(childParentList, thisGnChildren)
    
    sharedStats <- dplyr::bind_rows(sharedStats, thisGnParents)
  }
  
}


readr::write_tsv(
  x = sharedStats,
  file = file.path(outDir, "pangenome_phage_similarity.tab")
)


