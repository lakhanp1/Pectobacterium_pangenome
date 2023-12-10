suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/compare_hg_sets.R")
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
################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
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
  file = confs$analysis$prophages$preprocessing$files$mash_upgma
)

################################################################################
# # homology group annotations
# hgs <- unique(unlist(regionHgs$hgs))
# hgCog <- AnnotationDbi::select(
#   x = panOrgDb, keys = hgs,
#   columns = c("COG_id", "COG_description")
# ) %>% 
#   unique()
# 
# hgPfam <- AnnotationDbi::select(
#   x = panOrgDb, keys = hgs,
#   columns = c("pfam_description")
# ) %>% 
#   unique()


################################################################################
# prepare clusterjs JSON for a cluster/grp
grp <- clusterList$phage_grp_36

grpDnd <- ape::keep.tip(phy = mashTree, tip = grp$members) %>% 
  ape::as.hclust.phylo() %>% 
  as.dendrogram() %>% 
  dendextend::ladderize()

slots <- labels(grpDnd)[order.dendrogram(grpDnd)]

grpHgUnion <- regionList[grp$members] %>%
  purrr::map("hgs") %>% unlist() %>% unique()

clusterJsonDf <- NULL
linksJsonDf <- NULL
geneToGroup <- NULL

slotNum <- 0
# reg <- slots[1]

for (reg in slots) {
  
  # regionList[[reg]]
  genomeHgs <- AnnotationDbi::select(
    x = panOrgDb, keys = regionList[[reg]]$genomeId,
    columns = c(
      "GID", "genePos", "chr_id", "chr_name", "start", "end", "strand",
      "mRNA_key", "genePos", "mRNA_id", "COG_description", "pfam_description"
    ),
    keytype = "genomeId"
  )
  
  regHgs <- dplyr::left_join(
    x = tibble::tibble(hg = regionList[[reg]]$hgs, chr_name = regionList[[reg]]$chr),
    y = genomeHgs,
    by = c("hg" = "GID", "chr_name")
  ) %>%
    dplyr::arrange(chr_id, genePos)
  
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
    tidyr::replace_na(replace = list(COG = "-", PFAM = "-")) %>% 
    dplyr::mutate(
      strand = dplyr::if_else(strand == "-", -1, 1, 1),
      start = as.integer(start),
      end = as.integer(end)
    ) %>% 
    dplyr::rename(
      uid = mRNA_key, label = mRNA_id, chr = chr_id
    )
  
  regGenes <- dplyr::select(
    regHgs, uid, label, chr_name, chr, start, end, strand
  )
  
  # names should be a tibble type as it needs to be an {object} in JSON 
  # and not a [list of {objects}]
  regGenes$names <- dplyr::select(
    regHgs, COG, PFAM, genePos, genomeId, hg
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
      genes = c(uid, label, chr, start, end, strand, names)
    ) %>% 
    dplyr::select(
      uid = lociUid, name = lociName, start = lociStart, end = lociEnd,
      -chr_name, genes
    )
  
  clusterJsonDf <- dplyr::bind_rows(
    clusterJsonDf,
    tibble::tibble(
      uid = reg, name = reg, slot = slotNum,
      loci = list(lociDf)
    )
  )
  
  # add the links with previous cluster
  if(slotNum > 0){
    target <- clusterJsonDf$loci[slotNum] %>%
      purrr::map_dfr(.f = function(.x){.x$genes}) %>%
      tidyr::unnest(cols = names) %>%
      dplyr::rename_with(
        .fn = ~paste("target_", .x, sep = ""), .cols = -hg
      )

    query <- clusterJsonDf$loci[slotNum + 1] %>%
      purrr::map_dfr(.f = function(.x){.x$genes}) %>%
      tidyr::unnest(cols = names) %>%
      dplyr::rename_with(
        .fn = ~paste("query_", .x, sep = ""), .cols = -hg
      )
    
    ## IMPROVEMENT: replace this with end-to-end alignment of region HG sets
    hgMatches <- dplyr::left_join(
      x = target, y = query, by = "hg"
    ) %>% 
      tidyr::unite(col = uid, target_uid, query_uid, sep = "_link_", remove = FALSE)
    
    regLinks <- tibble::tibble(
      uid = hgMatches$uid,
      identity = 1, similarity = 1,
      target = dplyr::select(hgMatches, uid = target_uid, label = target_label),
      query = dplyr::select(hgMatches, uid = query_uid, label = query_label)
    )
    
    linksJsonDf <- dplyr::bind_rows(linksJsonDf, regLinks)

  }
  
  slotNum <- slotNum + 1
}

## make groups JSON
groupsJsonDf <- dplyr::group_by(geneToGroup, hg) %>% 
  dplyr::summarise(
    genes = list(genes),
    groupFreq = n(),
    hidden = "false"
  ) %>% 
  dplyr::arrange(desc(groupFreq)) %>% 
  dplyr::select(uid = hg, label = hg, everything())

groupsJsonDf <- dplyr::left_join(
  groupsJsonDf,
  y = tibble::tibble(
    groupFreq = 1:max(groupsJsonDf$groupFreq),
    colour = viridis::viridis(n = max(groupsJsonDf$groupFreq), option = "magma")
  ),
  by = "groupFreq"
)


mergedJsonDf <- tibble::tibble(
  clusters = list(clusterJsonDf),
  links = list(linksJsonDf),
  groups = list(groupsJsonDf)
)

jsonlite::toJSON(
  mergedJsonDf,
  dataframe = "rows", pretty = TRUE
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


