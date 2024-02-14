suppressMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressMessages(library(tidyverse))

# build HG PAV table for a few genomes of interest

rm(list = ls())

source("scripts/utils/homology_groups.R")

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

# genomes <- c("g_418", "g_162", "g_172")
genomes <- c("g_187", "g_399", "g_179", "g_191")

subPanHgs <- sub_pangenome_hgs(pandb = panOrgDb, genomes = genomes)

df <- AnnotationDbi::select(
  panOrgDb, keys = genomes,
  columns = c("genomeId", "GID", "COG_description", "pfam_description"),
  keytype = "genomeId"
) %>%
  dplyr::group_by(GID) %>% 
  dplyr::mutate(
    dplyr::across(
      .cols = c(COG_description, pfam_description),
      .fns = ~paste(unique(.x), collapse = "; ")
    ),
    value = 1
  ) %>% 
  dplyr::distinct() %>% 
  dplyr::ungroup() %>% 
  tidyr::pivot_wider(
    names_from = genomeId, values_fill = 0
  ) %>% 
  dplyr::left_join(y = subPanHgs, by = c("GID" = "hgId")) %>% 
  dplyr::arrange(desc(nGenomes)) %>% 
  dplyr::select(GID, nGenomes, class, !!!genomes, everything())

clipr::write_clip(df)
