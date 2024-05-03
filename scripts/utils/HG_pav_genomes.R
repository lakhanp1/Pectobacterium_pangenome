suppressMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressMessages(library(tidyverse))

# build HG PAV table for a few genomes of interest

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")

################################################################################

genomeOfInterest <- c("g_187", "g_404", "g_175")

# use phylogenetic tree nodes to categorize HGs into core, accessory and unique
coreClade <- "n59"
nodeExclude <- NULL

treeMethod <- "kmer_nj" # ani_upgma, kmer_upgma
panOrgDb <- org.Pectobacterium.spp.pan.eg.db

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

column_nodepath <- paste("nodepath.", treeMethod, sep = "")
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

subGenomes <- sampleInfo

if(!is.null(coreClade)){
  subGenomes <- dplyr::filter(
    sampleInfo,
    stringr::str_detect(
      string = !!sym(column_nodepath), pattern = paste(";", coreClade, ";", sep = "")
    )
  )
}

if(!is.null(nodeExclude)){
  subGenomes <- dplyr::filter(
    subGenomes,
    # stringr::str_detect(
    #   string = !!sym(column_nodepath), pattern = paste(";", coreClade, ";", sep = "")
    # ),
    stringr::str_detect(
      string = !!sym(column_nodepath), pattern = paste(";", nodeExclude, ";", sep = ""),
      negate = TRUE
    )
  )
}

################################################################################

# HG PAV for genomes of interest
genomeHgsPav <- sub_pangenome_hgs(pandb = panOrgDb, genomes = genomeOfInterest)

# HG position
hgPos <- AnnotationDbi::select(
  panOrgDb, keys = genomeOfInterest,
  columns = c("genomeId", "GID", "chr_id", "start", "end"),
  keytype = "genomeId"
) %>% 
  dplyr::mutate(
    dplyr::across(.cols = c(start, end), .fns = as.numeric),
  ) %>% 
  dplyr::group_by(genomeId) %>% 
  dplyr::arrange(chr_id, start) %>% 
  dplyr::mutate(hgPosition = 1:n()) %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(
    id_cols = c(GID),
    names_from = genomeId, values_from = hgPosition,
    values_fn = ~ min(.x)
  ) %>% 
  dplyr::select(hgId = GID, !!!genomeOfInterest) %>% 
  tidyr::replace_na(
    as.list(rep(-1, length(genomeOfInterest))) %>% purrr::set_names(genomeOfInterest)
  ) %>% 
  dplyr::arrange(!!!syms(genomeOfInterest))

# annotation data for genomes of interest
hgAn <- AnnotationDbi::select(
  panOrgDb, keys = genomeOfInterest,
  columns = c("genomeId", "GID", "COG_description", "pfam_description"),
  keytype = "genomeId"
) %>%
  dplyr::group_by(GID) %>% 
  dplyr::summarise(
    dplyr::across(
      .cols = c(COG_description, pfam_description),
      .fns = ~paste(unique(.x), collapse = "; ")
    )
  ) %>% 
  tidyr::unite(
    col = cog_pfam, COG_description, pfam_description, sep = "; ", remove = FALSE
  ) %>% 
  dplyr::mutate(
    broad_function = stringr::str_extract_all(
      string = cog_pfam,
      pattern = regex(
        "integrase|secretion.system|methyltransferase|restriction|integrase|transposase|toxin",
        ignore_case = TRUE
      )
    ),
    broad_function = purrr::map_chr(
      .x = broad_function,
      .f = function(x){
        if(length(x) > 0){
          return(paste(unique(tolower(x)), collapse = ";"))
        } else {
          return("")
        }
      }
    )
  ) %>% 
  dplyr::select(-cog_pfam) %>% 
  dplyr::relocate(COG_description, pfam_description, .after = broad_function)

# combine data
df <- dplyr::left_join(x = hgPos, y = genomeHgsPav, by = "hgId") %>% 
  dplyr::left_join(y = hgAn, by = c("hgId" = "GID"))

################################################################################
# add pangenome level and a sub-pangenome level core, accessory, unique annotation

# HG PAV for pangenome
panHgs <- sub_pangenome_hgs(
  pandb = panOrgDb, genomes = sampleInfo$genomeId
) %>% 
  dplyr::rename(pan_class = class, pan_nGenomes = nGenomes)

df <- dplyr::left_join(x = df, y = panHgs, by = "hgId")

if(!identical(sampleInfo, subGenomes)){
  # get HG PAV for genomes under a node 
  subPanHgs <- sub_pangenome_hgs(
    pandb = panOrgDb,
    genomes = union(genomeOfInterest, subGenomes$genomeId)
  ) %>% 
    dplyr::rename(subset_class = class, subset_nGenomes = nGenomes)
  
  df <- dplyr::left_join(x = df, y = subPanHgs, by = "hgId")

}

clipr::write_clip(df)

