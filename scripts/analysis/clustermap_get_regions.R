#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# 1) extract genomic coordinates between two homology groups of interest for genomes
# 2) prepare GFF3 file that include homology groups(as CDS) and associated PFAM and
# COG annotations

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
################################################################################
set.seed(124)

grpToView <- "ctv_pbr"

hgs <- c("hg_22427604", "hg_22427603")

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
treeMethod <- "kmer_upgma" # ani_upgma, kmer_upgma

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/", grpToView, sep = "")
outPrefix <- paste(outDir, "/", grpToView, sep = "")

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

cmJs <- jsonlite::read_json(
  path =  paste(outPrefix, ".json", sep = "")
)

genomes <- purrr::map_dfr(
  .x = cmJs$clusters,
  .f = function(x){
    aGene <- x$loci[[1]]$genes[[1]]
    list(
      prophage_id = x$name,
      genomeId = aGene$names$genomeId,
      chr = aGene$chr
    )
  }
) %>% 
  dplyr::mutate(inputOrder = 1:n())

hgPos <- AnnotationDbi::select(
  x = panOrgDb, keys = hgs,
  columns = c("genomeId", "chr_id", "chr_name", "start", "end", "strand")
) %>% 
  dplyr::mutate(
    dplyr::across(.cols = c(start, end), .fns = as.numeric)
  ) %>% 
  dplyr::right_join(y = genomes, by = c("chr_id" = "chr", "genomeId")) %>% 
  tibble::as_tibble() %>% 
  dplyr::add_count(genomeId, chr_id) %>% 
  dplyr::filter(n == 2) %>% 
  dplyr::arrange(genomeId, start) %>% 
  dplyr::mutate(
    regionStart = min(start),
    regionEnd = max(end),
    regionWidth = regionEnd - regionStart,
    .by = genomeId, chr_id,
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, genomeId, sampleId),
    by = "genomeId"
  ) %>% 
  dplyr::arrange(inputOrder)

regions <- dplyr::filter(hgPos, GID == hgs[1]) %>% 
  dplyr::mutate(
    region = paste(chr_name, ":", regionStart, "-", regionEnd, sep = "")
  ) %>% 
  dplyr::select(
    prophage_id, genomeId, sampleId, chr_name, regionStart, regionEnd, strand, region
  )

# generate GFF3 files for regions to later creation of GenBank files
for (i in 1:nrow(regions)) {
  regObj <- as.list(regions[i, ])
  
  # get gene annotation from panOrdDb
  regHgs <- region_homology_groups(
    pandb = panOrgDb, genome = regObj$genomeId,
    chr = regObj$chr_name, start = regObj$regionStart, end = regObj$regionEnd,
    cols = c(
      "GID", "genePos", "chr_id", "chr_name", "start", "end", "strand",
      "mRNA_key", "genePos", "mRNA_id", "COG_description", "pfam_description"
    )
  ) %>%
    dplyr::rename(ID = GID) %>% 
    dplyr::mutate(
      dplyr::across(.cols = c(start, end), .fns = as.numeric),
      locus_tag = ID,
      type = "CDS"
    ) %>% 
    dplyr::group_by(ID, start) %>% 
    dplyr::mutate(
      cog = paste(unique(COG_description), collapse = " "),
      pfam = paste(unique(pfam_description), collapse = " "),
      dplyr::across(
        .cols = c(cog, pfam),
        .fns = ~stringr::str_replace_all(string = .x, pattern = ",", replacement = " ::")
        )
    ) %>% 
    dplyr::slice(1L) %>% 
    dplyr::ungroup() %>% 
    dplyr::select(
      chr = chr_name, start, end, strand, type, ID, locus_tag, mRNA_id, mRNA_key, cog, pfam
    ) %>% 
    dplyr::arrange(start)
  
  # flip strand if needed
  if (regObj$strand == "-") {
    
    regHgs <- dplyr::mutate(
      regHgs,
      strand = dplyr::if_else(strand == "-", -1, 1, 1)
    ) %>% 
      invert_coordinates() %>% 
      dplyr::mutate(
        strand = dplyr::if_else(strand == -1, "-", "+", "+")
      ) %>% 
      dplyr::select(-starts_with("_"))
  }
  
  # store as GFF3
  regGr <- GenomicRanges::makeGRangesFromDataFrame(regHgs, keep.extra.columns = TRUE) %>% 
    GenomicRanges::shift(shift = -min(regHgs$start) + 1)
  
  rtracklayer::export.gff3(
    object = regGr,
    con = paste(outDir, "/region_fasta/", regObj$prophage_id, ".gff3", sep = "")
  )
  
}

dplyr::select(regions, prophage_id, genomeId, sampleId, region, strand) %>% 
  readr::write_tsv(
    file = paste(outPrefix, ".variable_regions.tab", sep = "")
  )

