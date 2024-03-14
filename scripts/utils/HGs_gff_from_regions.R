#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# 1) prepare GFF3 file that include homology groups(as CDS) and associated PFAM and
# COG annotations from a provided region
# 2) this GFF3 file will be used to generate a GenBank file which will be used
# in Mauve aligner

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser,
  opt_str = c("-r", "--regions"), type = "character", action = "store",
  help = "regions file (TSV)"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-g", "--genomes"), type = "character", action = "store",
  help = "optionally, a COMMA separated list of genome identifiers"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$regions))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# opts$regions <- "analysis/pangenome_v2/carotovoricin/ctv_hgt/ctv_region/hg_regions.tab"
# opts$genomes <- "g_145,g_194,g_429"

stopifnot(file.exists(opts$regions))
################################################################################

panOrgDb <- org.Pectobacterium.spp.pan.eg.db
outDir <- paste(dirname(opts$regions), "/gff3", sep = "")

if(!dir.exists(outDir)){ dir.create(outDir) }


hgRegions <- suppressMessages(readr::read_tsv(opts$regions))

if(!is.null(opts$genomes)){
  genomeSet <- stringr::str_split_1(string = opts$genomes, pattern = ",")
  hgRegions <- tibble::tibble(genomeId = intersect(genomeSet, hgRegions$genomeId)) %>% 
    dplyr::left_join(hgRegions, by = "genomeId")
}

stopifnot(
  c("chr_name", "start", "end", "strand", "genomeId", "regionId") %in% colnames(hgRegions),
  nrow(hgRegions) > 0
)

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
  
  regGr$phase <- 0
  
  rtracklayer::export.gff3(
    object = regGr,
    con = paste(outDir, "/", regObj$regionId, ".gff3", sep = "")
  )
}

