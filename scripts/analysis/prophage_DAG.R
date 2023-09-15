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

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

logger::log_threshold(DEBUG)
file_log <- paste(
  "logs/prophage_correction.", format(Sys.time(), "%Y%m%d.%H%M"), ".log",
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
  dplyr::select(prophage_id, prophage_length = length, chr)

# read prophage HGs stored locally
proHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$prophage_hg)
) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  ) %>%
  dplyr::rename(prophage_id = id) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, SpeciesName, genomeId, N50), by = "sampleId"
  ) %>%
  dplyr::left_join(y = prophageDf, by = "prophage_id") %>%
  # dplyr::filter(prophage_length >= 5000) %>%
  dplyr::arrange(desc(nHgs), desc(N50))

proHgL <- purrr::transpose(proHgs) %>%
  purrr::set_names(nm = purrr::map(., "prophage_id"))

sharedStats <- tibble::tibble()
childParentList <- list()

# create graph and add vertices to graph
phageGraph <- igraph::make_empty_graph() +
  igraph::vertices(proHgs$prophage_id, genomeId = proHgs$genomeId, nHgs = proHgs$nHgs)

# g_386.vir_2; g_438.vir_3

# sampleInfo <- dplyr::filter(
#   sampleInfo,
#   # genomeId %in% "g_378",
#   genomeId %in% c("g_93", "g_218", "g_3", "g_391", "g_438")
# )
# gn <- "g_278"
# ph <- "g_111.vir_3"
# ph <- "g_261.vir_5"

# iterating over genomes to detect parent prophages for each genomes's prophages
# this is important because a genome can have fragmented prophages because of
# breaks in assembly. Detecting parent for all prophages in the genome together
# will ensure that fragmented prophages in genome are assigned a same parent prophage
for (gn in sampleInfo$genomeId) {
  logger::log_info("processing phages in genome: ", gn)

  gnPhages <- dplyr::filter(proHgs, genomeId == !!gn) %>%
    dplyr::pull(prophage_id)

  logger::log_info("prophages in genome: ", paste(gnPhages, collapse = "; "))

  # purrr::keep(.x = proHgL, .p = function(x) x$genomeId == gn)

  gnChildPhages <- gnPhages
  bestParent <- list()

  # iterate over all phages to find the best parent phage
  for (ph in proHgs$prophage_id) {
    parentPh <- purrr::pluck(proHgL, !!!ph)

    if (parentPh$genomeId == gn) {
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

      # IMP: child -> parent asignment criteria
      if (nUniqHgs <= 2 && nSharedHgs > nUniqHgs) {
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


      childHgs <- purrr::map(proHgL[childPhages], "hgs") %>%
        unlist() %>%
        unique()
      
      # additional check for MULTIPLE childPhages --> parent relationship:
      # these childPhages will be merged to form a single prophage
      # 1) check that all phages are located on independent contigs or
      # at 5'/3' ends of large contigs
      contigCheck <- TRUE

      if (length(childPhages) > 1) {
        phageContigPos <- purrr::map_dbl(
          proHgL[childPhages],
          .f = ~ get_hg_sets_location(hgs = .x$hgs, genome = gn, chr = .x$chr, pandb = panOrgDb)
        )

        if (any((phageContigPos != 0) & (phageContigPos != 1))) {
          contigCheck <- FALSE
        }
      }

      # check for duplicate genome child-parent & parent-child circular links

      # collect data and comparison scores for the current childPhages --> parent
      thisParent <- list(
        parent = parentPh$prophage_id,
        parentGenome = parentPh$genomeId,
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
      thisParent$jaccardIndex <- thisParent$nSharedHgs / length(union(childHgs, parentPh$hgs))
      thisParent$contentDissimilarity <- 1 - mean(c(thisParent$perSharedParent, thisParent$perSharedChild))

      thisParent$relation <- dplyr::case_when(
        thisParent$jaccardIndex == 1 ~ "identical",
        thisParent$jaccardIndex >= 0.9 ~ "high_quality",
        thisParent$jaccardIndex >= 0.5 ~ "medium_quality",
        thisParent$jaccardIndex < 0.5 ~ "low_quality",
        contigCheck == FALSE ~ "NA",
        TRUE ~ "NA"
      )

      logger::log_debug("contigCheck: ", contigCheck)

      # compare with previous best parent
      if (contigCheck) {
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

        parentChildSetFound <- FALSE
        betterParentFound <- FALSE

        # iterate over parent set and replace the parent element combination
        ## if the current one is better
        for (i in 1:length(bestParent)) {
          # make sure correct set is being compared
          if (any(thisParent$children %in% bestParent[[i]]$children) && !betterParentFound) {
            parentChildSetFound <- TRUE

            oldParentPhage <- purrr::pluck(proHgL, !!!bestParent[[i]]$parent)
          
          logger::log_debug(
            "evaluating previous best parent: #", i, ": ", bestParent[[i]]$parent,
            " <-- ", bestParent[[i]]$childString,
            " [nHgs ", thisParent$nSharedHgs, "==", bestParent[[i]]$nSharedHgs,
            ", Len ", thisParent$childPhageLen, "==", bestParent[[i]]$childPhageLen,
            ", N50: ", parentPh$N50, " > ", oldParentPhage$N50, "]"
          ) 
          
          logger::log_debug(
            "parentChildSetFound: ", parentChildSetFound,
            ", betterParentFound: ", betterParentFound
          )

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
              logger::log_info(
                "Found better parent: sharedHgs = ", compSharedHgs,
                ", length = ", compLen, " :", isBetterParent
              )
              logger::log_info(
                "nHgs ", thisParent$nSharedHgs, "==", bestParent[[i]]$nSharedHgs,
                ", Len ", thisParent$childPhageLen, "==", bestParent[[i]]$childPhageLen,
                ", N50: ", parentPh$N50, " > ", oldParentPhage$N50
              )
              logger::log_info(
                "Previous parent: #", i, " : ",
                paste(bestParent[[i]]$children, collapse = ";"), " --parent--> ",
                bestParent[[i]]$parent
              )
              logger::log_info(
                "New parent: #", i, " : ",
                paste(thisParent$children, collapse = ";"), " --parent--> ",
                thisParent$parent
              )

              betterParentFound <- TRUE
              bestParent[[i]] <- thisParent

              logger::log_info("#parents: ", length(bestParent))
            }
          } else if (betterParentFound) {
            # if better parent was found in one of the previous iteration,
            # ensure the other parent-child relationship that involve any
            # of the child in this iteration is set to NA

            if (any(thisParent$children %in% bestParent[[i]]$children)) {
              logger::log_info(
                "Removing parent #", i,
                " because betterParentFound =", betterParentFound, " in previous iteration: ",
                paste(bestParent[[i]]$children, collapse = ";"), " --parent--> ",
                bestParent[[i]]$parent
              )

              bestParent[[i]] <- NULL

              logger::log_info("#parents: ", length(bestParent))
            }
          }
        }

        # parent exists but not for this child phages, add the current phage to bestParentList
        if (!parentChildSetFound) {
          bestParent <- append(bestParent, list(thisParent))

          logger::log_info(
            "New parent: #", length(bestParent), " because ",
            "parentChildSetFound = ", parentChildSetFound, ": ",
            paste(thisParent$children, collapse = ";"), " --parent--> ",
            thisParent$parent
          )
          logger::log_info("#parents: ", length(bestParent))
        }
      }
    }
  }

  thisGnParents <- NULL

  if (length(bestParent) > 0) {
    thisGnParents <- purrr::map_dfr(
      .x = bestParent, .f = `[`,
      c(
        "childString", "childGenome", "childPhageLen", "nHgChild",
        "parent", "parentGenome", "nHgParent", "nSharedHgs", 
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
    thisGnParents <- purrr::map_dfr(
      .x = proHgL[noParents], .f = `[`,
      c("prophage_id", "genomeId", "nHgs")
    ) %>%
      dplyr::rename(child = prophage_id, childGenome = genomeId, nHgChild = nHgs) %>%
      dplyr::bind_rows(thisGnParents)
  }

  # finally merge with sharedStats global data
  sharedStats <- dplyr::bind_rows(sharedStats, thisGnParents)
}


readr::write_tsv(
  x = sharedStats,
  file = confs$analysis$prophages$files$phage_similarity
)


