#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(logger))
suppressPackageStartupMessages(library(configr))

# for all prophage fragments in a genome, find a best parent prophage
# if two or more prophage fragments in genome can be linked to a single larger
# prophage, use this to merge prophages

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
prophageLenCutoff <- confs$analysis$prophages$cutoff_length
maxSyntenyGap <- 2
minSyntenyChain <- 5

outDir <- confs$analysis$prophages$preprocessing$dir

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

logger::log_threshold(DEBUG)
file_log <- paste(
  "logs/prophage_merge.", format(Sys.time(), "%Y%m%d.%H%M"), ".log",
  sep = ""
)
logger::log_appender(
  appender = logger::appender_file(file = file_log)
)

################################################################################
# sorting by N50 to avoid the circular links. In case of duplicate prophages
# across multiple genome assemblies, prophage from assembly with highest N50 is
# always assigned as parent to the prophages from remaining assembly. Prophage
# from the highest quality assembly is root in this case
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) %>%
  dplyr::arrange(N50)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

################################################################################

prophageDf <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>%
  dplyr::select(prophage_id, prophage_length, chr)

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$raw_prophage_hg)
) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName, genomeId, N50), by = "sampleId"
  ) %>%
  dplyr::left_join(y = prophageDf, by = "prophage_id") %>%
  dplyr::arrange(desc(nHgs), desc(N50))

proHgL <- purrr::transpose(proHgs) %>%
  purrr::set_names(nm = purrr::map(., "prophage_id"))

sharedStats <- tibble::tibble()
childParentList <- list()

# create graph and add vertices to graph
phageGraph <- igraph::make_empty_graph() +
  igraph::vertices(proHgs$prophage_id, genomeId = proHgs$genomeId, nHgs = proHgs$nHgs)

# sampleInfo <- dplyr::filter(
#   sampleInfo,
#   # genomeId %in% "g_378",
#   genomeId %in% c("g_93", "g_218", "g_3", "g_391", "g_438")
# )
# gn <- "g_278"
# gn <- "g_128"

# iterating over genomes to detect parent prophages for each genomes's prophages
# this is important because a genome can have fragmented prophages because of
# breaks in assembly. Detecting parent for all prophages in the genome together
# will ensure that fragmented prophages in genome are assigned a same parent prophage
for (gn in sampleInfo$genomeId) {
  logger::log_info("processing phages in genome: ", gn)
  
  gnPhages <- dplyr::filter(proHgs, genomeId == !!gn) %>%
    dplyr::pull(prophage_id)
  
  logger::log_info("prophages in genome: ", paste(gnPhages, collapse = "; "))
  
  gnChildPhages <- gnPhages
  bestParent <- list()
  
  # iterate over all phages to find the best parent phage
  for (ph in proHgs$prophage_id) {
    parentPh <- purrr::pluck(proHgL, !!!ph)
    
    # do not compare with very small prophages or prophages from same genome
    if (parentPh$genomeId == gn || parentPh$prophage_length < prophageLenCutoff) {
      next
    }
    
    childPhages <- NULL
    
    # for each prophage in genome, check for overlap with prophage list
    # this will help to identify a union of prophages in genome that overlap
    # with one parent prophage
    for (gp in gnChildPhages) {
      gpL <- purrr::pluck(proHgL, !!!gp)
      
      nUniqHgs <- length(setdiff(gpL$hgs, parentPh$hgs))
      nSharedHgs <- length(intersect(gpL$hgs, parentPh$hgs))
      
      logger::log_debug(
        gp, " --parent--> ", ph, " || nUniqHgs: ",
        nUniqHgs, "; nSharedHgs: ", nSharedHgs
      )
      
      # IMP: child -> parent assignment criteria
      if (nUniqHgs <= 5 && nSharedHgs > nUniqHgs) {
        childPhages <- append(childPhages, gpL$prophage_id)
        logger::log_debug(gp, " --parent--> ", ph, ": significant match found")
      }
    }
    
    # evaluate if this parent is the best, if not, update the parent
    if (!is.null(childPhages)) {
      logger::log_info(
        paste(childPhages, collapse = "; "), " --parent--> ", ph,
        " : evaulating relationship"
      )
      
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
      
      if (all(isLinkedPreviously == 1)) {
        logger::log_warn(
          paste(childPhages, collapse = "; "), " --parent--> ", parentPh$prophage_id,
          " : Circular link avoided"
        )
        
        next
      }
      
      
      # get the shared syntenic homology groups between child and parent
      syntenicOverlap <- purrr::map(
        .x = proHgL[childPhages],
        .f = ~ syntenic_hg_overlap(
          ref = parentPh$hgs, qur = .x$hgs,
          minChainLen = min(length(.x$hgs), minSyntenyChain),
          maxGapLen = maxSyntenyGap
        )
      )
      
      # check if there is very large gaps between phage fragments aligned to parent phage
      # need
      # to
      # develop
      # this
      # improvement
      
      # only consider child phages those share syntenic HGs
      # remove child phages that did not share any syntenic HGs
      syntenicOverlap <- purrr::compact(syntenicOverlap)
      childPhages <- names(syntenicOverlap)
      
      # skip child phages if no syntenic shared HGs
      syntenicShared <- purrr::map(syntenicOverlap, "lcs") %>%
        unlist() %>%
        unique()
      
      if (length(syntenicShared) == 0) {
        next
      }
      
      childPhageLen <- purrr::map_dbl(proHgL[childPhages], "prophage_length") %>% sum()

      
      # shared HGs between child and parent
      childHgs <- purrr::map(proHgL[childPhages], "hgs") %>%
        unlist() %>%
        unique()
      
      sharedHgs <- intersect(childHgs, parentPh$hgs)
      
      # order of fragmented phages on the parent phage
      orderedChildPhages <- purrr::map_dfr(
        .x = syntenicOverlap,
        .f = function(x){
          list(
            pos = list(x$pos$s1),
            start = x$pos$s1[1],
            end = tail(x$pos$s1, n = 1)
          )
        },
        .id = "child"
      ) %>% 
        dplyr::arrange(start)
      
      childPos <- dplyr::select(orderedChildPhages, child, start, end) %>% 
        tidyr::unite(col = "fragmentAlnPos", start, end, sep = ":") %>% 
        dplyr::summarise(
          dplyr::across(
            .cols = everything(),
            .fns = ~paste(., collapse = ";")
          )
        )

      
      # additional check for MULTIPLE childPhages --> parent relationship:
      # these childPhages will be merged to form a single prophage
      # 1) check that all phages are located on independent contigs or
      # at 5'/3' ends of large contigs
      fragmentedPhageContigs <- TRUE
      
      if (length(childPhages) > 1) {
        phageContigPos <- purrr::map_dbl(
          proHgL[childPhages],
          .f = ~ get_hg_sets_location(
            hgs = .x$hgs, genome = gn, chr = .x$chr, pandb = panOrgDb
          )
        )
        
        logger::log_debug(
          "Fragmented prophages with positions: ",
          stringi::stri_join(names(phageContigPos), phageContigPos, collapse = "; ", sep = ": ")
        )
        
        if (any((phageContigPos != 0) & (phageContigPos != 1))) {
          fragmentedPhageContigs <- FALSE
        }
      }
      
      logger::log_debug("fragmentedPhageContigs: ", fragmentedPhageContigs)
      
      # collect data and comparison scores for the current childPhages --> parent
      thisParent <- list(
        parent = parentPh$prophage_id,
        parentGenome = parentPh$genomeId,
        nHgParent = parentPh$nHgs,
        children = childPhages,
        childGenome = gn,
        childString = childPos$child,
        fragmentAlnPos = childPos$fragmentAlnPos,
        nChildFragments = length(childPhages),
        childHgs = childHgs,
        nHgChild = length(childHgs),
        nSharedHgs = length(sharedHgs),
        nSyntenicSharedHgs = length(syntenicShared),
        nUniqParent = length(setdiff(parentPh$hgs, syntenicShared)),
        nUniqChild = length(setdiff(childHgs, syntenicShared)),
        childPhageLen = childPhageLen
      )
      
      thisParent$perSharedParent <- thisParent$nSyntenicSharedHgs / thisParent$nHgParent
      thisParent$perSharedChild <- thisParent$nSyntenicSharedHgs / thisParent$nHgChild
      thisParent$jaccardIndex <- thisParent$nSyntenicSharedHgs / length(union(childHgs, parentPh$hgs))
      thisParent$contentDissimilarity <- 1 - mean(c(thisParent$perSharedParent, thisParent$perSharedChild))
      
      thisParent$relation <- dplyr::case_when(
        thisParent$jaccardIndex == 1 ~ "identical",
        thisParent$jaccardIndex >= 0.9 | thisParent$contentDissimilarity < 0.1 ~ "high_quality",
        thisParent$jaccardIndex >= 0.5 ~ "medium_quality",
        thisParent$jaccardIndex < 0.5 ~ "low_quality",
        fragmentedPhageContigs == FALSE ~ "NA",
        TRUE ~ "NA"
      )
      
      # compare with previous best parent and update if the current one is a better match
      if (fragmentedPhageContigs) {
        # if no parent found until now: make current phage as parent
        if (length(bestParent) == 0) {
          bestParent[1] <- list(thisParent)
          logger::log_info(
            thisParent$childString, " --parent--> ", thisParent$parent,
            ": First parent"
          )
          logger::log_info("#parents: ", length(bestParent))
          
          next
        }
        
        parentExistForThisChild <- FALSE
        betterParentFound <- FALSE
        
        betterParentIndex <- c()
        matchedExistingParents <- c()
        
        logger::log_info("#parents: ", length(bestParent))
        
        # if there are parent detected previously:
        # iterate over parent set and replace the parent element combination
        # if the current one is better
        for (i in 1:length(bestParent)) {
          oldParentPhage <- purrr::pluck(proHgL, !!!bestParent[[i]]$parent)
          parentStatus <- NULL
          
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
          
          # evaluate the current parent with previous parent set
          if (any(thisParent$children %in% bestParent[[i]]$children)) {
            parentExistForThisChild <- TRUE
            
            matchedExistingParents <- append(matchedExistingParents, i)
            
            if (isBetterParent != "no") {
              betterParentFound <- TRUE
              betterParentIndex <- append(betterParentIndex, i)
            }
          }
          
          parentStatus <- paste("isBetterParent:", isBetterParent)
          
          logger::log_debug(
            "Evaluated previous best parent: #", i, ": ", bestParent[[i]]$parent,
            " <--parent-- ", bestParent[[i]]$childString,
            " [nHgs ", thisParent$nSharedHgs, "==", bestParent[[i]]$nSharedHgs,
            ", Len ", thisParent$childPhageLen, "==", bestParent[[i]]$childPhageLen,
            ", N50: ", parentPh$N50, " > ", oldParentPhage$N50, "]",
            "|| action: ", parentStatus
          )
        }
        
        # warn if multiple parents have child in common
        childToManyParents <- FALSE
        
        if (length(matchedExistingParents) > 1) {
          logger::log_info(
            "matchedExistingParents = #", paste(matchedExistingParents, collapse = ";"), ": ",
            purrr::map_chr(bestParent[matchedExistingParents], "parent") %>%
              paste(collapse = ";")
          )
          
          existingChildSet <- purrr::map(bestParent[matchedExistingParents], "children") %>%
            unlist() %>%
            table()
          
          childToManyParents <- any(existingChildSet > 1)
          
          if (childToManyParents) {
            purrr::map(
              bestParent[matchedExistingParents],
              ~ logger::log_warn(
                "Detected multiple parents that share child: ",
                .x$parent, "  <--parent--  ", .x$childString
              )
            )
          }
        }
        
        
        # adjust the parent list if a better parent was identified
        if (betterParentFound) {
          logger::log_info(
            "betterParentIndex = #", paste(betterParentIndex, collapse = ";"), ": ",
            purrr::map_chr(bestParent[betterParentIndex], "parent") %>%
              paste(collapse = ";")
          )
          
          if (setequal(betterParentIndex, matchedExistingParents)) {
            # remove all parents that have child from thisParent$child
            bestParent[betterParentIndex] <- NULL
            bestParent <- append(bestParent, list(thisParent))
            
            logger::log_info(
              "Found better parent compared to the existing: #: ", thisParent$parent,
              " <--parent-- ", paste(thisParent$children, collapse = ";")
            )
          } else {
            logger::log_info(
              "Found a better parent but cannot replace because ",
              "betterParentIndex != matchedExistingParents"
            )
          }
        }
        
        
        # parent/s exists but not for this child phages, add the current phage to bestParentList
        if (!parentExistForThisChild) {
          bestParent <- append(bestParent, list(thisParent))
          
          logger::log_info(
            "Found a new parent-child set: #: ", thisParent$parent,
            " <--parent-- ", paste(thisParent$children, collapse = ";")
          )
        }
        
        logger::log_info("#parents: ", length(bestParent))
      }
    }
  }
  
  thisGnParents <- NULL
  
  if (length(bestParent) > 0) {
    thisGnParents <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c(
        "childString", "nChildFragments", "childGenome", "childPhageLen",
        "fragmentAlnPos", "nHgChild",
        "parent", "parentGenome", "nHgParent", "nSharedHgs", "nSyntenicSharedHgs",
        "jaccardIndex", "contentDissimilarity",
        "perSharedParent", "perSharedChild", "relation"
      )
    ) %>%
      dplyr::rename(child = childString)
    
    edgesDf <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c(
        "children", "parent", "nSharedHgs", "jaccardIndex", "contentDissimilarity",
        "perSharedParent", "perSharedChild", "relation"
      )
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
    distinctPhages <- purrr::map_dfr(
      .x = proHgL[noParents], .f = `[`,
      c("prophage_id", "genomeId", "nHgs", "prophage_length")
    ) %>%
      dplyr::rename(
        child = prophage_id, childGenome = genomeId,
        nHgChild = nHgs, childPhageLen = prophage_length
      ) %>% 
      dplyr::mutate(nChildFragments = 1)
    
    thisGnParents <- dplyr::bind_rows(thisGnParents, distinctPhages)
  }
  
  if(!is.null(thisGnParents)){
    
    thisGnParents <- dplyr::mutate(
      thisGnParents, keyCol = 1:n()
    )
    
    # give a new identifier to the merged prophage groups
    if (any(thisGnParents$nChildFragments > 1)){
      
      mergedPhages <- dplyr::filter(
        thisGnParents, nChildFragments > 1
      ) %>% 
        dplyr::arrange(desc(childPhageLen)) %>% 
        dplyr::mutate(
          prophage_id = paste(childGenome, ".merged_vir_", 1:n(), sep = ""),
          .before = child
        ) %>% 
        dplyr::select(prophage_id, keyCol)
      
      thisGnParents <- dplyr::full_join(
        thisGnParents, mergedPhages, by = "keyCol"
      )
      
    } else {
      thisGnParents %<>%  dplyr::mutate(
        prophage_id = child
      )
    }
    
    # filter flag for very smaller child phages (eg. <5kb)
    thisGnParents <- dplyr::mutate(
      thisGnParents,
      prophage_id = dplyr::if_else(is.na(prophage_id), child, prophage_id),
      filtered = dplyr::if_else(
        condition = childPhageLen < prophageLenCutoff, true = 1, false = 0
      )
    ) %>% 
      dplyr::select(prophage_id, filtered, everything(), -keyCol) %>% 
      dplyr::rename(
        fragments = child, nFragments = nChildFragments, genomeId = childGenome,
        nHg = nHgChild, prophage_length = childPhageLen
      )
    
    # finally merge with sharedStats global data
    sharedStats <- dplyr::bind_rows(sharedStats, thisGnParents)
  }
  
}


readr::write_tsv(
  x = sharedStats,
  file = confs$analysis$prophages$preprocessing$files$consolidated_phages
)


