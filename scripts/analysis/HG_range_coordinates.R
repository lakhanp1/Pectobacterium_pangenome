#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))

# 1) extract genomic coordinates between two homology groups of interest for genomes
# 2) prepare GFF3 file that include homology groups(as CDS) and associated PFAM and
# COG annotations
# 3) this GFF3 file will be used to generate a GenBank file which will be used
# in Mauve aligner

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
################################################################################
set.seed(124)

grpName <- "ctv_conserved"

excludeHgs <- FALSE
makeGff <- FALSE

# hgs <- c("hg_22427643", "hg_22427599")      # ctv + flanking 3 genes
hgs <- c("hg_22427640", "hg_22427604")    # ctv conserved loci
# hgs <- c("hg_22427604", "hg_22427603")    # ctv tail loci

# genomes <- c(
#   "g_145", "g_194", "g_429", "g_442", "g_421", "g_426",
#   "g_150", "g_447", "g_434",
#   "g_57", "g_53", "g_221", "g_125",
#   "g_158", "g_446", "g_66", "g_222", "g_296", "g_442", "g_8", "g_38", "g_273", "g_259",
#   "g_305", "g_378", "g_428", "g_248", "g_449", "g_54", "g_116", "g_423", "g_375", "g_381"
# )

genomes <- NULL

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- paste(confs$analysis$prophages$dir, "/cluster_viz/ctv_hgt/", grpName, sep = "")

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

if(is.null(genomes) | length(genomes) == 0){
  genomes <- sampleInfo$genomeId
}

hgPos <- AnnotationDbi::select(
  x = panOrgDb, keys = hgs,
  columns = c(
    "GID", "mRNA_id",
    "genomeId", "chr", "chr_id", "chr_name", "start", "end", "strand"
  )
) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(
    dplyr::across(.cols = c(start, end), .fns = as.numeric),
    length = end - start + 1
  ) %>% 
  dplyr::filter(genomeId %in% !!genomes) %>% 
  dplyr::add_count(genomeId, chr_id) %>%
  dplyr::filter(n == 2)

hgGrl <- GenomicRanges::makeGRangesListFromDataFrame(
  df = hgPos, split.field = "genomeId", keep.extra.columns = TRUE
)

hgRange <- unlist(range(hgGrl, ignore.strand = TRUE))
hgRegionGr <- hgRange

if(excludeHgs){
  hgRegionGr <- GenomicRanges::psetdiff(
    x = hgRange, y = hgGrl, ignore.strand = TRUE
  ) %>% 
    unlist()
}

hgRegions <- as.data.frame(hgRegionGr) %>% 
  tibble::rownames_to_column(var = "genomeId") %>% 
  dplyr::select(genomeId, chr = seqnames, start, end) %>% 
  tibble::as_tibble() %>% 
  dplyr::left_join(
    y = dplyr::filter(hgPos, GID == !!hgs[1]) %>% 
      dplyr::select(genomeId, chr, chr_name, strand),
    by = c("genomeId", "chr")
  ) %>%
    dplyr::left_join(
      y = dplyr::select(sampleInfo, genomeId, sampleId, SpeciesName),
      by = "genomeId"
    ) %>% 
  dplyr::mutate(
    region = paste(chr_name, ":", start, "-", end, sep = ""),
    SpeciesName = stringr::str_replace_all(
      string = SpeciesName, pattern = "\\.?\\s+", replacement = "_"
    ),
    regionId = paste(genomeId, ".", SpeciesName, sep = "")
  ) %>%
  dplyr::select(
    regionId, genomeId, sampleId, chr_name, start, end, strand, region
  )


################################################################################
# generate GFF3 files for hgRegions to later creation of GenBank files
if(makeGff){
  for (i in 1:nrow(hgRegions)) {
    regObj <- as.list(hgRegions[i, ])
    
    # get gene annotation from panOrdDb
    regHgs <- region_homology_groups(
      pandb = panOrgDb, genome = regObj$genomeId,
      chr = regObj$chr_name, start = regObj$start, end = regObj$end,
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
          .fns = ~ stringr::str_replace_all(string = .x, pattern = ",", replacement = " ::")
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
      con = paste(outDir, "/", regObj$regionId, ".gff3", sep = "")
    )
  }
}

dplyr::select(hgRegions, regionId, genomeId, sampleId, region, strand) %>%
  readr::write_tsv(
    file = paste(outDir, "/", "hg_regions.tab", sep = "")
  )


