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

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- confs$analysis$prophages$dir

orgDb <- org.Pectobacterium.spp.pan.eg.db

logger::log_threshold(WARN)
file_log <- "logs/prophage_correction.log"
unlink(file_log)
logger::log_appender(
  appender = logger::appender_file(file = file_log)
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

# create graph and add vertices to graph
phageGraph <- igraph::make_empty_graph() + 
  igraph::vertices(proHgs$prophage_id, Genome = proHgs$Genome, nHgs = proHgs$nHgs)

# g_386.vir_2; g_438.vir_3

# sampleInfo <- dplyr::filter(
#   sampleInfo,
#   # Genome %in% "378",
#   Genome %in% c("68", "87", "94", "90", "59", "86", "30", "93")
# )
# gn <- "378"
# ph <- "g_111.vir_3"
# ph <- "g_261.vir_5"

# iterating over genomes to detect parent prophages for each genomes's prophages
# this is important because a genome can have fragmented prophages because of
# breaks in assembly. Detecting parent for all prophages in the genome together
# will ensure that fragmented prophages in genome are assigned a same parent prophage
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
    # this will help to identify a union of prophages in genome that overlap
    # with one parent prophage
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
      
      logger::log_debug('parent phage: ', ph, ' | evaulating children: ',
                        paste(childPhages, collapse = "; "))
      
      # prevent circular links: 
      # was there a previous link from parentPh-(has_parent)->childPhages
      # no need to make childPhages-(has_parent)->parentPh link
      # to the current childPhages
      isLinkedPreviously <- purrr::map_dbl(
        .x = childPhages,
        .f = ~ igraph::edge_connectivity(
          graph = phageGraph, source = parentPh$prophage_id, target = .x
        )
      )
      
      if(all(isLinkedPreviously == 1)){
        logger::log_warn('Circular link avoided between ', parentPh$prophage_id,
                         ' and ', paste(childPhages, collapse = "; "))
        
        break
      }
      
      
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
            
            logger::log_debug(
              "parent: #", i, ": ", bestParent[[i]]$parent,
              "nHgs ", thisParent$nSharedHgs, "==", bestParent[[i]]$nSharedHgs,
              ", Len ", thisParent$childPhageLen, "==", bestParent[[i]]$childPhageLen,
              ", N50: ", parentPh$N50, " > ", oldParentPhage$N50
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
  
  thisGnParents <- NULL
  
  if (length(bestParent) > 0) {
    
    thisGnParents <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c("childString", "childGenome", "nHgChild", "parent", "parentGenome", "nHgParent",
        "nSharedHgs", "perSharedParent", "perSharedChild", "relation")
    ) %>% 
      dplyr::rename(child = childString)
    
    edgesDf <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c("children", "parent", "nSharedHgs", "perSharedParent",
        "perSharedChild", "relation")
    )
    
    phageEdges <- purrr::map2(
      .x = edgesDf$children, .y = edgesDf$parent, .f = c
    ) %>% 
      unlist()
    
    # add edges with attributes: direction is child_phage--(has_parent)-->parent_phage
    phageGraph <- phageGraph + 
      igraph::edges(
        phageEdges,
        nSharedHgs = edgesDf$nSharedHgs, perSharedChild = edgesDf$perSharedChild,
        perSharedParent = edgesDf$perSharedParent, relation = edgesDf$relation
      )
    
    # list for circular child-parent-child linking handling
    thisGnChildren <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c("children", "parent")
    ) %>% 
      tibble::deframe() %>% 
      as.list()
    
    childParentList <- append(childParentList, thisGnChildren)
    
  }
  
  # for remaining phages where parent was not found
  noParents <- setdiff(gnPhages, unlist(purrr::map(bestParent, "children")))
  if (length(noParents) > 0) {
    thisGnParents <- purrr::map_dfr(
      .x = proHgL[noParents], .f = `[`,
      c("prophage_id", "Genome", "nHgs")
    ) %>% 
      dplyr::rename(child = prophage_id, childGenome = Genome, nHgChild = nHgs) %>% 
      dplyr::bind_rows(thisGnParents)
  }
  
  # finally merge with sharedStats global data
  sharedStats <- dplyr::bind_rows(sharedStats, thisGnParents)
  
}


readr::write_tsv(
  x = sharedStats,
  file = confs$analysis$prophages$files$prophage_similarity 
)


