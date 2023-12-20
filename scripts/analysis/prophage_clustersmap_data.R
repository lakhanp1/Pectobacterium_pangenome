suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
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

grpToView <- "phage_grp_1"
subSample <- TRUE
################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
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
  dplyr::left_join(y = rawRegions, by = "prophage_id") %>% 
  dplyr::left_join(y = regionHgs, by = "prophage_id")

################################################################################
regionList <- purrr::transpose(regionClusters) %>% 
  purrr::set_names(nm = purrr::map(., "prophage_id"))

clusterList <- dplyr::group_by(regionClusters, phage_grp) %>% 
  dplyr::group_map(
    .f = ~{
      list(
        phage_grp = .x$phage_grp[1],
        members = .x$prophage_id,
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
  mixed = c("CI repressor", "hydrolase", "holin", "portal protein",
            "VRR-NUC domain", "other"),
  replication = c("DNA polymerase", "reverse transcriptase"),
  tail = c("tail fibre", "tail tube", "tail", "tail collar", "baseplate",
           "tail sheath", "tail tape"),
  capsid = c("capsid", "head"),
  defense = c("colicin", "toxin-antitoxin"),
  unknown = c("unknown", "uncharacterized conserved")
)

hgFuncColors <- purrr::map2(
  .x = funcTypes,
  .y = c(viridis::viridis(n = length(funcTypes) -1, option = "turbo"), "grey"),
  .f = function(x, y){
    tibble::tibble(
      function_category = x,
      function_color_hex = y,
      function_color = paste(
        "rgb(",
        paste(as.vector(col2rgb(y)), collapse = ", "),
        ")", sep = "")
    )
  }
) %>% 
  purrr::list_rbind(names_to = "broad_function")

################################################################################
# prepare clusterjs JSON for a cluster/grp
grp <- clusterList[[grpToView]]

# frequency for all HGs in the prophages in current cluster
grpHgFreq <- regionList[grp$members] %>%
  purrr::map("hgs") %>% unlist() %>% table() %>% 
  tibble::enframe(name = "hgId", value = "freq") %>% 
  dplyr::mutate(freq = as.numeric(freq)) %>% 
  dplyr::left_join(phageHgTypes, by = "hgId") %>% 
  dplyr::left_join(hgFuncColors, by = "function_category") %>% 
  dplyr::arrange(desc(freq))

# original MASH distance tree
grpDnd <- ape::keep.tip(phy = mashTree, tip = grp$members) %>% 
  ape::as.hclust.phylo() %>% 
  as.dendrogram() %>% 
  dendextend::ladderize()

plot(rev(grpDnd), horiz = TRUE)

# HG PAV clustering
hgPavMat <- purrr::map(
  .x = regionList[labels(grpDnd)],
  .f = function(x){
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

hgPavDnd <- hclust(d = dist(hgPavMat)) %>% as.dendrogram() %>% 
  dendextend::ladderize()

plot(rev(hgPavDnd), horiz = TRUE)
# order.dendrogram(hgPavDnd) == sort(hgPavDnd, type = "nodes") %>% order.dendrogram() %>% rev()


if(subSample){
  
  grpSubset <- dendextend::cutree(tree = hgPavDnd, h = 2) %>% 
    tibble::enframe(name = "prophage_id", value = "cut") %>% 
    dplyr::add_count(cut, name = "count") %>% 
    dplyr::slice_sample(n = 1, by = cut)
  
  subHgPavDnd <- dendextend::prune(
    dend = hgPavDnd,
    leaves = setdiff(labels(hgPavDnd), grpSubset$prophage_id)
  ) %>% 
    dendextend::ladderize()
  
  grpSubset <- dplyr::left_join(
    x = grpSubset,
    y = tibble::tibble(prophage_id = labels(subHgPavDnd), order = 1:nleaves(subHgPavDnd)),
    by = "prophage_id"
  ) %>% 
    dplyr::arrange(order)
  
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
  slots <- labels(grpDnd)
}

hgStrand <- NULL
clusterJsonDf <- NULL
linksJsonDf <- NULL
geneToGroup <- NULL

slotNum <- 0
# reg <- slots[1]

for (reg in slots) {
  
  regObj <- regionList[[reg]]
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
      dplyr::across(.cols = c(start, end, genePos), .fns = as.integer)
    ) %>% 
    dplyr::rename(
      uid = mRNA_key, label = mRNA_id, chr = chr_id
    ) %>%
    dplyr::arrange(start)
  
  # change the orientation of gene strand for better visualization
  thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
    as.list()
  
  if(!is.null(hgStrand)){
    
    for (h in grpHgFreq$hgId) {
      
      if(!is.null(thisRegHgStrand[[h]]) & !is.null(hgStrand[[h]])){
        if((thisRegHgStrand[[h]] != hgStrand[[h]])){
          
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
    dplyr::select(regHgs, hg, genes = uid)
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
groupsJsonDf <- dplyr::group_by(geneToGroup, hg) %>% 
  dplyr::summarise(
    genes = list(genes),
    groupFreq = n(),
    hidden = FALSE
  ) %>% 
  dplyr::arrange(desc(groupFreq)) %>% 
  dplyr::select(uid = hg, label = hg, everything()) %>% 
  dplyr::left_join(y = phageHgTypes, by = c("uid" = "hgId"))

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
  groupsJsonDf, hgFreqColors, by = "groupFreq"
) %>% 
  dplyr::left_join(hgFuncColors, by = "function_category") %>% 
  dplyr::rename(colour = function_color) %>% 
  dplyr::arrange(broad_function) %>% 
  dplyr::arrange(desc(groupFreq))

## make links JSON
linksJsonDf <- purrr::map2_dfr(
  .x = groupsJsonDf$uid,
  .y = groupsJsonDf$genes,
  .f = function(h, g){
    if(length(g) > 1){
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
  path = paste(outDir, "/cluster_viz/clustermap_config.json", sep = "")
)

# unbox the json for correct format
unbox_list_json <- function(x){
  purrr::map(
    .x = x,
    .f = function(y){
      if(is.list(y)){
        y <- unbox_list_json(y)
      } else{
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
  path = paste(outDir, "/cluster_viz/", grpToView, ".json", sep = "")
)

################################################################################

## progressive alignment test
# masterHgs <- proHgL[[slotOrder[1]]]$hgs
# 
# i <- 3
# seq1 <- masterHgs
# seq2 <- proHgL[[slotOrder[i]]]$hgs
# 
# lcs <- syntenic_hg_overlap(ref = seq1, qur = seq2)
# 
# if (lcs$strand$s2 == "-") {
#   seq2 <- rev(seq2)
# }


