#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# combine DefenseFinder, PHROG, COG, PFAM annotations for prophage homology groups

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

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

prophagePool <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$prophage_pool)
)

phage_genes <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$raw_prophage_genes)
)

################################################################################

## import DefenseFinder annotations
defense_finder <- suppressMessages(
  readr::read_tsv(
    confs$data$defense_finder$files$defense_systems,
    col_names = c(
      "sampleId", "sys_id", "type", "subtype", "sys_beg", "sys_end",
      "protein_in_syst", "genes_count", "name_of_profiles_in_sys"
    )
  )
)

defense_systems <- dplyr::left_join(
  dplyr::select(sampleInfo, sampleId, genomeId, SpeciesName),
  defense_finder,
  by = "sampleId"
)

defense_genes <- suppressMessages(
  readr::read_tsv(
    confs$data$defense_finder$files$defense_genes,
    col_names = c(
      "sampleId", "replicon", "hit_id", "gene_name", "hit_pos", "model_fqn",
      "sys_id", "sys_loci", "locus_num", "sys_wholeness", "sys_score", "sys_occ",
      "hit_gene_ref", "hit_status", "hit_seq_len", "hit_i_eval", "hit_score",
      "hit_profile_cov", "hit_seq_cov", "hit_begin_match", "hit_end_match",
      "counterpart", "used_in"
    )
  )
) %>%
  dplyr::select(sampleId, sys_id, mRNA_id = hit_id, defense_gene = gene_name) %>%
  dplyr::full_join(defense_systems, by = c("sampleId", "sys_id")) %>%
  dplyr::select(genomeId, sys_id, type, mRNA_id, defense_gene)

phage_defense <- dplyr::left_join(
  x = dplyr::select(phage_genes, genomeId, prophage_id, mRNA_id),
  y = defense_genes,
  by = c("genomeId", "mRNA_id")
) %>%
  dplyr::filter(!is.na(sys_id)) %>%
  dplyr::group_by(sys_id) %>%
  dplyr::summarise(
    prophage_id = paste(unique(prophage_id), collapse = ";"),
    prophage_defense_system = TRUE
  )

defense_systems <- dplyr::left_join(
  defense_systems, phage_defense,
  by = c("sys_id")
) %>%
  tidyr::replace_na(list(prophage_defense_system = FALSE))

readr::write_tsv(
  x = defense_systems,
  file = confs$analysis$defense_systems$files$data
)

################################################################################
## import PHROG annotations
phrogs_meta <- suppressMessages(
  readr::read_tsv(confs$data$phrog$files$phrog_metadata)
) %>%
  dplyr::mutate(phrog = paste("phrog", phrog, sep = "_")) %>%
  dplyr::select(phrog, phrog_annotation = annot, phrog_category = category)

## non-redundant protein cluster mappings
nr_clusters <- suppressMessages(
  readr::read_tsv(
    confs$data$pangenomes$pectobacterium.v2$proteins$files$nr_protein_clusters,
    col_names = c("nr_cluster", "member")
  )
)

## non-redundant pan-proteome phrog annotation
nr_phrogs <- suppressMessages(
  readr::read_tsv(
    file = confs$data$phrog$files$annotation,
    col_names = c(
      "phrog", "seqId", "alnScore", "seqIdentity", "eVal",
      "qStart", "qEnd", "qLen", "tStart", "tEnd", "tLen")
  )
)

## best phrog hit for NR-pan-proteome
nr_phrog_best_hit <- dplyr::group_by(nr_phrogs, seqId) %>%
  dplyr::arrange(eVal, .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

## map NR-pan-proteome phrog annotations to each member of the cluster
phrog_best_hit <- dplyr::left_join(
  x = nr_clusters, y = nr_phrog_best_hit,
  by = c("nr_cluster" = "seqId")
) %>%
  dplyr::rename(seqId = member) %>%
  tidyr::separate_wider_regex(
    cols = seqId, cols_remove = FALSE,
    patterns = c(genomeId = "^\\d+", "_", mRNA_id = ".*$")
  ) %>%
  dplyr::mutate(genomeId = paste("g_", genomeId, sep = ""))

phage_phrogs <- dplyr::left_join(
  x = phage_genes, y = phrog_best_hit,
  by = c("genomeId", "mRNA_id")
) %>%
  dplyr::select(
    prophage_id, nr_cluster, seqId, genomeId, mRNA_id, chr_name, start, end,
    strand, hgId, phrog, seqIdentity, eVal, qLen, tLen
  ) %>%
  dplyr::left_join(y = phrogs_meta, by = "phrog") %>%
  dplyr::add_count(hgId, name = "hg_count") %>%
  dplyr::arrange(desc(hg_count), hgId, phrog_category)

################################################################################
## consolidate DefenseFinder, PHROG, COG, PFAM annotations at the homology group level

## DefenseFinder
hg_defense <- dplyr::left_join(
  x = dplyr::select(phage_genes, genomeId, prophage_id, mRNA_id, hgId),
  y = defense_genes,
  by = c("genomeId", "mRNA_id")
) %>%
  dplyr::filter(!is.na(sys_id)) %>%
  dplyr::select(hgId, defense_system = type, defense_gene) %>%
  dplyr::distinct() %>%
  dplyr::group_by(hgId) %>%
  dplyr::reframe(
    dplyr::across(
      .cols = starts_with("defense"),
      .fns = ~ paste(.x, collapse = ";")
    )
  )

## PHROGs
hg_phrog <- dplyr::filter(phage_phrogs, !is.na(phrog)) %>%
  dplyr::select(hgId, starts_with("phrog_")) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(phrog_annotation)) %>%
  dplyr::group_by(hgId) %>%
  dplyr::reframe(
    dplyr::across(
      .cols = starts_with("phrog"),
      .fns = ~ paste(.x, collapse = ";")
    )
  ) %>%
  dplyr::mutate(annotation_db = "PHROG")

hg_cog <- suppressMessages(
  AnnotationDbi::select(
    x = panOrgDb, keys = unique(phage_genes$hgId),
    columns = c("GID", "COG_id", "COG_description")
  )
) %>%
  dplyr::filter(COG_id != "None") %>%
  dplyr::distinct() %>%
  dplyr::group_by(GID) %>%
  dplyr::reframe(
    dplyr::across(
      .cols = starts_with("COG_"),
      .fns = ~ paste(.x, collapse = ";")
    )
  )

hg_pfam <- suppressMessages(
  AnnotationDbi::select(
    x = panOrgDb, keys = unique(phage_genes$hgId),
    columns = c("GID", "pfam_id", "pfam_description")
  )
) %>%
  dplyr::filter(!is.na(pfam_id)) %>%
  dplyr::distinct() %>%
  dplyr::group_by(GID) %>%
  dplyr::reframe(
    dplyr::across(
      .cols = starts_with("pfam_"),
      .fns = ~ paste(.x, collapse = ";")
    )
  )

phageHgDf <- suppressMessages(
  AnnotationDbi::select(
    x = panOrgDb, keys = unique(phage_genes$hgId), columns = c("GID", "class")
  )
) %>%
  dplyr::distinct() %>%
  dplyr::rename(hgId = GID) %>%
  # dplyr::left_join(hg_phrog, by = c("hgId")) %>%
  dplyr::left_join(hg_cog, by = c("hgId" = "GID")) %>%
  dplyr::left_join(hg_pfam, by = c("hgId" = "GID")) %>%
  dplyr::left_join(hg_defense, by = c("hgId")) %>%
  dplyr::mutate(
    annotation = NA_character_,
    annotation_category = NA_character_
  ) %>%
  dplyr::select(
    hgId, annotation, annotation_category,
    starts_with("defense_"), starts_with("COG"), starts_with("pfam"), class
  )

readr::write_tsv(
  x = phageHgDf,
  file = confs$analysis$prophages$preprocessing$files$hg_annotations
)

