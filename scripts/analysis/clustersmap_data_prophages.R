suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
source("scripts/utils/heatmap_utils.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

grpToView <- "phage_grp_1"
subSample <- TRUE
cutHeight <- 1.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: host phylogeny, prophage HG PAV, prophage MASH,
# completeness score
clusterOrder <- "host" # host, hg_pav, cluster_mash

# a vector of prophage identifiers that will be included in clustermap plot
appendPhages <- c()

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list()

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
treeMethod <- "kmer_upgma" # ani_upgma, kmer_upgma

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

rawTree <- import_tree(
  file = confs$analysis$phylogeny[[treeMethod]]$files$tree,
  phylo = TRUE
)

hostDnd <- as.dendrogram(ape::as.hclust.phylo(rawTree)) %>%
  dendextend::ladderize()

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny[[treeMethod]]$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>% dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

rawRegions <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>%
  dplyr::select(prophage_id, chr, start, end, topology)

# read prophage HGs stored locally
regionHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$raw_prophage_hg)
) %>%
  dplyr::select(prophage_id, hgs) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  )

regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
) %>%
  # handle fragmented regions in the future
  # dplyr::mutate(fragments = stringr::str_split(fragments, ";")) %>%
  # tidyr::unnest(fragments) %>%
  dplyr::left_join(y = rawRegions, by = "prophage_id") %>%
  dplyr::left_join(y = regionHgs, by = "prophage_id")

################################################################################
regionList <- purrr::transpose(regionClusters) %>%
  purrr::set_names(nm = purrr::map(., "prophage_id"))

regionList <- append(regionList, customRegions)

clusterList <- dplyr::group_by(regionClusters, phage_grp) %>%
  dplyr::group_map(
    .f = ~ {
      list(
        phage_grp = .x$phage_grp[1],
        members = .x$prophage_id[.x$nFragments == 1],
        # save fragmented regions separately for now.
        fragmented = .x$prophage_id[.x$nFragments > 1],
        group_size = nrow(.x)
      )
    },
    .keep = TRUE
  ) %>%
  purrr::set_names(nm = purrr::map(., "phage_grp"))

mashTree <- ape::read.tree(
  file = confs$analysis$prophages$preprocessing$files$mash_hclust
)

phageHgTypes <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$hg_broad_functions)
)

funcTypes <- list(
  DNA = c("nuclease", "helicase", "terminase", "DNA function"),
  integrase = "integrase",
  mixed = c(
    "CI repressor", "hydrolase", "holin", "portal protein",
    "VRR-NUC domain", "other"
  ),
  replication = c("DNA polymerase", "reverse transcriptase"),
  tail = c(
    "tail fibre", "tail tube", "tail", "tail collar", "baseplate",
    "tail sheath", "tail tape"
  ),
  capsid = c("capsid", "head"),
  defense = c("colicin", "toxin-antitoxin"),
  unknown = c("unknown", "uncharacterized conserved"),
  flanking = c("flanking")
)

hgFuncColors <- purrr::map2(
  .x = funcTypes,
  .y = c(viridis::viridis(n = length(funcTypes) - 2, option = "turbo"), "grey50", "white"),
  .f = function(x, y) {
    tibble::tibble(
      function_category = x,
      function_color_hex = y,
      function_color = paste(
        "rgb(",
        paste(as.vector(col2rgb(y)), collapse = ", "),
        ")",
        sep = ""
      )
    )
  }
) %>%
  purrr::list_rbind(names_to = "broad_function")

################################################################################
# prepare clusterjs JSON for a cluster/grp
grp <- clusterList[[grpToView]]

# # optionally, a custom region list can be provided to generate the plot
# grp <- list(
#   phage_grp = grpToView,
#   members = c( ),
#   fragmented = character(),
#   group_size = integer()
# )

# grp <- list(
#   phage_grp = grpToView,
#   members = dplyr::filter(
#     regionClusters,
#     SpeciesName == "P. brasiliense", nFragments == 1, phage_grp == "phage_grp_1",
#     !prophage_id %in% c("g_408.vir_3", "g_403.vir_3", "g_399.vir_3")
#   ) %>%
#     dplyr::pull(prophage_id)
# )

outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/", grp$phage_grp, sep = "")
outPrefix <- paste(outDir, "/", grp$phage_grp, sep = "")

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# frequency for all HGs in the prophages in current cluster
grpHgFreq <- regionList[grp$members] %>%
  purrr::map("hgs") %>%
  unlist() %>%
  table() %>%
  tibble::enframe(name = "hgId", value = "freq") %>%
  dplyr::mutate(freq = as.numeric(freq)) %>%
  dplyr::left_join(phageHgTypes, by = "hgId") %>%
  dplyr::left_join(hgFuncColors, by = "function_category") %>%
  dplyr::arrange(desc(freq))

# original MASH distance tree
grpMashDnd <- ape::keep.tip(phy = mashTree, tip = grp$members) %>%
  ape::as.hclust.phylo() %>%
  as.dendrogram() %>%
  dendextend::ladderize()

plot(rev(grpMashDnd), horiz = TRUE)

# HG PAV clustering
hgPavMat <- purrr::map(
  .x = regionList[labels(grpMashDnd)],
  .f = function(x) {
    table(x$hgs) %>%
      tibble::enframe(name = "hgs", value = "count") %>%
      dplyr::mutate(count = as.numeric(count))
  }
) %>%
  purrr::list_rbind(names_to = "prophage_id") %>%
  tidyr::pivot_wider(
    id_cols = prophage_id, names_from = hgs,
    values_from = count, values_fill = 0
  ) %>%
  tibble::column_to_rownames(var = "prophage_id") %>%
  as.matrix()

hgPavDnd <- hclust(d = dist(hgPavMat)) %>%
  as.dendrogram() %>%
  dendextend::ladderize()

plot(rev(hgPavDnd), horiz = TRUE)
# order.dendrogram(hgPavDnd) == sort(hgPavDnd, type = "nodes") %>% order.dendrogram() %>% rev()


# # plot HG PAV matrix with host taxonomy
# memberHostData <- tibble::tibble(genomeId = grp$members) %>% 
#   dplyr::left_join(
#     y = dplyr::select(regionClusters, genomeId = prophage_id, SpeciesName),
#     by = "genomeId"
#   )
# 
# spHt <- species_key_heatmap(
#   genomes = rownames(hgPavMat), speciesInfo = memberHostData
# )
# 
# spHt@heatmap_param$width <- unit(10, "cm")
# spHt@row_dend_param$width <- unit(6, "cm")
# 
# ht <- homology_group_heatmap(
#   mat = hgPavMat, phy = hgPavDnd,
#   width = unit(10, "cm"),
#   name = "HG_PAV", column_title = "HG PAV",
#   use_raster = TRUE, raster_quality = 3
# )
# 
# # ht@column_dend_param$cluster <- FALSE
# # ht@column_names_param$show <- FALSE
# 
# htList <- spHt + ht
# 
# ComplexHeatmap::draw(
#   object = htList,
#   main_heatmap = "HG_PAV",
#   row_dend_side = "left",
#   merge_legends = TRUE,
#   heatmap_legend_side = "bottom"
# 
# )


grpMemberData <- dplyr::left_join(
  x = tibble::tibble(prophage_id = c(grp$members, appendPhages)),
  y = dplyr::select(regionClusters, prophage_id, phage_grp, genomeId, SpeciesName,
                    isolation_source, geo_loc_country),
  by = "prophage_id"
) %>% 
  dplyr::left_join(
    y = tibble::tibble(
      genomeId = labels(hostDnd),
      host_order = 1:dendextend::nleaves(hostDnd)
    ),
    by = "genomeId"
  ) %>%
  dplyr::left_join(
    y = tibble::tibble(
      prophage_id = labels(grpMashDnd), mash_order = 1:dendextend::nleaves(grpMashDnd)
    ),
    by = "prophage_id"
  ) %>%
  dplyr::left_join(
    y = tibble::tibble(
      prophage_id = labels(hgPavDnd), pav_order = 1:dendextend::nleaves(hgPavDnd)
    ),
    by = "prophage_id"
  )


if (subSample) {
  pavDndCut <- dendextend::cutree(tree = hgPavDnd, h = cutHeight) %>%
    tibble::enframe(name = "prophage_id", value = "cut") %>%
    dplyr::add_count(cut, name = "count")
  
  grpSubset <- dplyr::slice_sample(pavDndCut, n = 1, by = cut)
  
  if (length(appendPhages) > 0) {
    grpSubset <- dplyr::bind_rows(
      grpSubset,
      tibble::tibble(prophage_id = appendPhages)
    ) %>%
      dplyr::distinct(prophage_id, .keep_all = TRUE)
  }
  
  grpMemberData <- dplyr::left_join(
    x = grpMemberData, y = pavDndCut, by = "prophage_id"
  ) %>%
    dplyr::left_join(
      y = tibble::tibble(prophage_id = grpSubset$prophage_id, shown = 1),
      by = "prophage_id"
    ) %>% 
    tidyr::replace_na(replace = list(shown = 0))
  
  
  subHgPavDnd <- dendextend::prune(
    dend = hgPavDnd,
    leaves = setdiff(labels(hgPavDnd), grpSubset$prophage_id)
  ) %>%
    dendextend::ladderize()
  
  subMashDnd <- dendextend::prune(
    dend = grpMashDnd,
    leaves = setdiff(labels(grpMashDnd), grpSubset$prophage_id)
  ) %>%
    dendextend::ladderize()
  
  grpSubset <- dplyr::left_join(
    x = grpSubset,
    y = tibble::tibble(
      prophage_id = labels(subHgPavDnd),
      order = 1:nleaves(subHgPavDnd)
    ),
    by = "prophage_id"
  ) %>%
    dplyr::arrange(order)
  
  hgPavDnd <- dendextend::set(
    hgPavDnd,
    what = "labels_col",
    value = as.numeric(labels(hgPavDnd) %in% grpSubset$prophage_id) + 1
  )
  
  grpMashDnd <- dendextend::set(
    grpMashDnd,
    what = "labels_col",
    value = as.numeric(labels(grpMashDnd) %in% grpSubset$prophage_id) + 1
  )
  
  plot(subHgPavDnd, horiz = TRUE)
  
  # sort(subHgPavDnd, type = "nodes") %>% order.dendrogram() %>% rev() == order.dendrogram(subHgPavDnd)
  # sort(subHgPavDnd, type = "nodes") %>% labels() %>% rev() == labels(subHgPavDnd)
  # labels(subHgPavDnd)
  # labels(subHgPavDnd)[order.dendrogram(subHgPavDnd)]
  # sort(subHgPavDnd, type = "nodes") %>% labels() %>% rev()
  # order.dendrogram(subHgPavDnd)
  # sort(subHgPavDnd, type = "nodes") %>% order.dendrogram()
  
  slots <- rev(grpSubset$prophage_id)
} else {
  slots <- labels(grpMashDnd)
  
  grpMemberData <- dplyr::mutate(grpMemberData, shown = 1)
}


readr::write_tsv(
  dplyr::arrange(grpMemberData, pav_order),
  file = paste(outPrefix, ".shown_clusters.txt", sep = "")
)

# show only specific groups if subsampled
grpMemberData <- dplyr::filter(grpMemberData, shown == 1)

viewClusters <- dplyr::case_when(
  clusterOrder == "host" ~
    grpMemberData$prophage_id[order(grpMemberData$host_order)],
  clusterOrder == "hg_pav" ~
    grpMemberData$prophage_id[order(grpMemberData$pav_order)],
  clusterOrder == "cluster_mash" ~
    grpMemberData$prophage_id[order(grpMemberData$mash_order)],
  TRUE ~ slots
)

viewClusters <- append(viewClusters, names(customRegions))

hgStrand <- NULL
clusterJsonDf <- NULL
linksJsonDf <- NULL
geneToGroup <- NULL
flankingHgs <- character(0)
slotNum <- 0
# reg <- viewClusters[1]

for (reg in viewClusters) {
  regObj <- regionList[[reg]]
  
  if (addFlankingRegions) {
    regObj$oldStart <- regObj$start
    regObj$oldEnd <- regObj$end
    
    if (!is.na(regObj$start)) {
      regObj$start <- pmax(regObj$oldStart - flankingRegion, 1)
    }
    
    if (!is.na(regObj$end)) {
      regObj$end <- regObj$end + flankingRegion
    }
  }
  
  regHgs <- region_homology_groups(
    pandb = panOrgDb, genome = regObj$genomeId,
    chr = regObj$chr, start = regObj$start, end = regObj$end,
    cols = c(
      "GID", "genePos", "chr_id", "chr_name", "start", "end", "strand",
      "mRNA_key", "genePos", "mRNA_id", "COG_description", "pfam_description"
    )
  ) %>%
    dplyr::rename(hg = GID)
  
  # summarize COG and PFAM annotations
  cog <- dplyr::select(regHgs, mRNA_key, COG_description) %>%
    dplyr::distinct() %>%
    dplyr::group_by(mRNA_key) %>%
    dplyr::summarise(COG = paste(COG_description, collapse = ";"))
  
  pfam <- dplyr::select(regHgs, mRNA_key, pfam_description) %>%
    dplyr::filter(!is.na(pfam_description)) %>%
    dplyr::distinct() %>%
    dplyr::group_by(mRNA_key) %>%
    dplyr::summarise(PFAM = paste(pfam_description, collapse = ";"))
  
  regHgs <- dplyr::select(regHgs, -COG_description, -pfam_description) %>%
    dplyr::distinct() %>%
    dplyr::left_join(y = cog, by = "mRNA_key") %>%
    dplyr::left_join(y = pfam, by = "mRNA_key") %>%
    dplyr::left_join(y = phageHgTypes, by = c("hg" = "hgId")) %>%
    tidyr::replace_na(replace = list(COG = "-", PFAM = "-")) %>%
    dplyr::mutate(
      strand = dplyr::if_else(strand == "-", -1, 1, 1),
      dplyr::across(.cols = c(start, end, genePos), .fns = as.integer),
      flanking = 0
    ) %>%
    dplyr::rename(
      uid = mRNA_key, label = mRNA_id, chr = chr_id
    ) %>%
    dplyr::arrange(start)
  
  # set function_category = flanking to the flanking genes added
  if (addFlankingRegions) {
    flankingGenes <- which(regHgs$start < regObj$oldStart | regHgs$end > regObj$oldEnd)
    
    regHgs$function_category[flankingGenes] <- "flanking"
    # regHgs$hg[flankingGenes] <- paste(
    #   regHgs$hg[flankingGenes], "_flanking", sep = ""
    # )
    
    regHgs$flanking[flankingGenes] <- 1
  }
  
  # change the orientation of gene strand for better visualization
  thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
    as.list()
  
  if (!is.null(hgStrand)) {
    for (h in grpHgFreq$hgId) {
      if (!is.null(thisRegHgStrand[[h]]) & !is.null(hgStrand[[h]])) {
        if ((thisRegHgStrand[[h]] != hgStrand[[h]])) {
          regHgs <- invert_coordinates(regHgs)
          
          # update this region strands backup record after changing orientation
          thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
            as.list(thisRegHgStrand)
        }
        break
      }
    }
  } else {
    for (h in grpHgFreq$hgId) {
      if (!is.null(thisRegHgStrand[[h]])) {
        if (thisRegHgStrand[[h]] == -1) {
          regHgs <- invert_coordinates(regHgs)
          
          # update this region strands backup record after changing orientation
          thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
            as.list(thisRegHgStrand)
        }
        
        break
      }
    }
  }
  
  hgStrand <- thisRegHgStrand
  
  # regHgs <- regHgs[1:5, ]
  
  regGenes <- dplyr::select(
    regHgs, uid, label, chr_name, chr, start, end, strand
  )
  
  # names should be a tibble type as it needs to be an {object} in JSON
  # and not a [list of {objects}]
  regGenes$names <- dplyr::select(
    regHgs, COG, PFAM, function_category, hg, genePos, genomeId
  )
  
  ## save groups information for making groups JSON
  geneToGroup <- dplyr::bind_rows(
    geneToGroup,
    dplyr::select(regHgs, hg, genes = uid, flanking)
  )
  
  # for now there is only one loci in each cluster
  # in future, there can be multiple loci when a region is combination of
  # multiple fragments
  lociDf <- dplyr::mutate(
    regGenes,
    lociUid = chr, lociName = chr_name,
    lociStart = min(start),
    lociEnd = max(end)
  ) %>%
    tidyr::nest(
      genes = c(uid, label, names, chr, start, end, strand)
    ) %>%
    dplyr::select(
      uid = lociUid, name = lociName, start = lociStart, end = lociEnd,
      -chr_name, genes
    )
  
  clusterJsonDf <- dplyr::bind_rows(
    clusterJsonDf,
    tibble::tibble(
      uid = stringr::str_replace_all(reg, pattern = "\\.", replacement = "_"),
      name = reg, slot = slotNum,
      regionId = reg,
      loci = list(lociDf)
    )
  )
  
  slotNum <- slotNum + 1
}

## make groups JSON
groupsJsonDf <- dplyr::mutate(
  geneToGroup,
  uid = dplyr::if_else(flanking == 1, paste(hg, "_flanking", sep = ""), hg)
) %>%
  dplyr::add_count(hg, name = "groupFreq") %>%
  dplyr::summarise(
    genes = list(genes),
    hidden = FALSE,
    .by = c(uid, hg, flanking, groupFreq)
  ) %>%
  dplyr::arrange(desc(groupFreq)) %>%
  dplyr::select(uid, label = hg, everything()) %>%
  dplyr::left_join(y = phageHgTypes, by = c("uid" = "hgId")) %>%
  dplyr::mutate(
    function_category = dplyr::if_else(
      flanking == 1, "flanking", function_category
    )
  )

# color HGs by their frequency
hgFreqColors <- tibble::tibble(
  groupFreq = 1:max(groupsJsonDf$groupFreq),
  freq_color_hex = viridis::viridis(
    n = max(groupsJsonDf$groupFreq), option = "magma", direction = -1
  )
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(freq_color = paste(as.vector(col2rgb(freq_color_hex)), collapse = ", ")) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(freq_color = paste("rgb(", freq_color, ")", sep = "")) %>%
  dplyr::select(-freq_color_hex)

groupsJsonDf <- dplyr::left_join(
  groupsJsonDf, hgFreqColors,
  by = "groupFreq"
) %>%
  dplyr::left_join(hgFuncColors, by = "function_category") %>%
  dplyr::rename(colour = function_color) %>%
  dplyr::arrange(broad_function) %>%
  dplyr::arrange(desc(groupFreq))

## make links JSON
linksJsonDf <- purrr::map2_dfr(
  .x = groupsJsonDf$uid,
  .y = groupsJsonDf$genes,
  .f = function(h, g) {
    if (length(g) > 1) {
      linkCombs <- combn(x = g, m = 2) %>%
        t() %>%
        as.data.frame.matrix() %>%
        dplyr::rename(
          query = V1, target = V2
        ) %>%
        tidyr::unite(col = uid, query, target, sep = "_link_", remove = F) %>%
        dplyr::mutate(
          identity = 1,
          similarity = 1,
          hg = h
        )
      
      regLinks <- dplyr::select(
        linkCombs, uid, identity, similarity
      ) %>%
        dplyr::mutate(
          target = dplyr::select(linkCombs, uid = target),
          query = dplyr::select(linkCombs, uid = query)
        )
      
      return(regLinks)
    }
  }
)

mergedJsonList <- list(
  clusters = clusterJsonDf,
  links = linksJsonDf,
  groups = groupsJsonDf
)

confJson <- jsonlite::read_json(
  path = paste(confs$analysis$prophages$dir, "/cluster_viz/clustermap_config.json", sep = "")
)

# unbox the json for correct format
unbox_list_json <- function(x) {
  purrr::map(
    .x = x,
    .f = function(y) {
      if (is.list(y)) {
        y <- unbox_list_json(y)
      } else {
        return(jsonlite::unbox(y))
      }
    }
  )
}

confJson <- unbox_list_json(confJson)

jsonlite::write_json(
  # x = list(config = confJson$config, data = mergedJsonList),
  x = mergedJsonList,
  # pretty = 2,
  path = paste(outPrefix, ".json", sep = "")
)

################################################################################
# prepare homology group PAV matrix from pan.db

markGenomeColors <- list(
  prophages = list(
    genomes = purrr::map_chr(.x = regionList[grp$members], .f = "genomeId") %>%
      unname(),
    color = "#D81B60"
  )
)

if (length(grp$fragmented) > 0) {
  markGenomeColors[["fragmented"]] <- list(
    genomes = purrr::map_chr(.x = regionList[grp$fragmented], .f = "genomeId") %>%
      unname(),
    color = "#FFC107"
  )
}

htSpecies <- species_key_heatmap(
  genomes = rawTree$tip.label, speciesInfo = sampleInfo,
  markGenomes = markGenomeColors
)

htSpecies@heatmap_param$width <- unit(12, "cm")

hgMat <- homology_groups_mat(
  pandb = panOrgDb, type = "pav", groups = grpHgFreq$hgId
)

hgMat <- hgMat[rawTree$tip.label, ]

funcTypeColors <- dplyr::select(
  hgFuncColors,
  at = broad_function, fill = function_color_hex
) %>%
  dplyr::distinct()

hgsInPan <- AnnotationDbi::select(
  x = panOrgDb, keys = grpHgFreq$hgId,
  columns = c(
    "genome", "genomeId", "chr_id", "chr_name",
    "start", "end", "strand", "mRNA_id", "genePos"
  )
) %>%
  dplyr::mutate(
    dplyr::across(.cols = c(start, end, genome, genePos), .fns = as.numeric)
  ) %>%
  dplyr::mutate(length = end - start) %>%
  dplyr::arrange(genome, chr_id, start)

readr::write_tsv(
  hgsInPan,
  file = paste(outPrefix, ".hgs_pangenome_data.txt", sep = "")
)

hgLengths <- split(x = hgsInPan$length, f = hgsInPan$GID)
hgLengths <- hgLengths[grpHgFreq$hgId]

# HG frequency barplot annotation
anFreq <- ComplexHeatmap::HeatmapAnnotation(
  length = ComplexHeatmap::anno_boxplot(
    x = hgLengths, which = "column"
  ),
  freq = ComplexHeatmap::anno_barplot(
    x = grpHgFreq$freq, which = "column"
  ),
  func = grpHgFreq$broad_function,
  col = list(
    func = deframe(funcTypeColors)
  ),
  gp = gpar(col = "black"),
  annotation_height = c(4, 4, 1), height = unit(3, "cm"),
  gap = unit(1, "mm"),
  border = c(length = TRUE, freq = TRUE, func = TRUE)
)

# HG PAV heatmap using PAV data for clustering
ht <- homology_group_heatmap(
  mat = hgMat, phy = TRUE,
  # width = unit(10, "cm"),
  name = "hgs",
  top_annotation = anFreq,
  use_raster = TRUE, raster_quality = 3
)

ht@column_dend_param$cluster <- FALSE
ht@column_names_param$show <- TRUE

htList <- htSpecies + ht

pdf(
  file = paste(outPrefix, ".hg_heatmap.pdf", sep = ""),
  width = 16, height = 10
)

ComplexHeatmap::draw(
  object = htList,
  column_title = paste(
    nrow(grpHgFreq), "homology groups for prophages in cluster:", grpToView
  ),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  main_heatmap = "hgs",
  row_dend_side = "left",
  merge_legends = TRUE,
  heatmap_legend_side = "right"
)

# use the phylogeny for clustering genomes
ComplexHeatmap::draw(
  object = htList,
  column_title = paste(
    nrow(grpHgFreq), "homology groups for prophages in cluster:", grpToView
  ),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  main_heatmap = "hgs",
  row_dend_side = "left",
  merge_legends = TRUE,
  cluster_rows = as.dendrogram_ordered.phylo(
    phy = rawTree, sourceOrder = rawTree$tip.label
  ),
  heatmap_legend_side = "right"
)

dendextend::tanglegram(
  dend1 = rev(grpMashDnd), main_left = "MASH",
  dend2 = rev(hgPavDnd), main_right = "HG PAV",
  lwd = 2, edge.lwd = 2, lab.cex = 1, margin_inner = 5, margin_outer = 2,
  common_subtrees_color_branches = TRUE,
  highlight_distinct_edges = FALSE,
  highlight_branches_lwd = FALSE
)

if (subSample) {
  dendextend::tanglegram(
    dend1 = rev(subMashDnd), main_left = "MASH",
    dend2 = rev(subHgPavDnd), main_right = "HG PAV",
    lwd = 2, edge.lwd = 2, lab.cex = 1, margin_inner = 5, margin_outer = 2,
    common_subtrees_color_branches = TRUE,
    highlight_distinct_edges = FALSE,
    highlight_branches_lwd = FALSE
  )
}

dev.off()


