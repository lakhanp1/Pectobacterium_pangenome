#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(optparse))

# 1) extract genomic coordinates between two homology groups of interest for genomes

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()


parser <- optparse::add_option(
  parser,
  opt_str = c("--hgs"), type = "character", action = "store",
  help = "a COMMA separated list of homology groups"
)

parser <- optparse::add_option(
  parser, default = NULL,
  opt_str = c("--genomes"), type = "character", action = "store",
  help = "a COMMA separated list of genome ids"
)

parser <- optparse::add_option(
  parser, default = FALSE,
  opt_str = c("--inner_region"), type = "logical", action = "store_true",
  help = "LOGICAL: TRUE: use the inner boundry of the region, othewise outer"
)


parser <- optparse::add_option(
  parser, default = FALSE,
  opt_str = c("--overlapping"), type = "logical", action = "store_true",
  help = "LOGICAL: If TRUE, include the overlapping gene's HG"
)

parser <- optparse::add_option(
  parser, default = FALSE,
  opt_str = c("--haplotypes"), type = "logical", action = "store_true",
  help = "LOGICAL: TRUE: optionally create haplotypes based on the homology
  groups signatures for all extracted regions"
)

parser <- optparse::add_option(
  parser, default = ".",
  opt_str = c("-d", "--dir"), type = "character", action = "store",
  help = "output dir"
)


parser <- optparse::add_option(
  parser, default = NULL,
  opt_str = c("-o", "--out"), type = "character", action = "store",
  help = "Output file name"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$hgs), is.null(opts$dir))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# opts$hgs <- "hg_22427604,hg_22427603"
# opts$dir <- "analysis/pangenome_v2/carotovoricin/ctv_tail"
# opts$inner_region <- TRUE
# opts$haplotypes <- TRUE
################################################################################

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

hgs <- stringr::str_split_1(string = opts$hgs, pattern = ",")

if(!dir.exists(opts$dir)){
  dir.create(opts$dir)
}

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

if(is.null(opts$genomes)){
  genomes <- sampleInfo$genomeId
} else{
  genomes <- stringr::str_split_1(string = opts$genomes, pattern = ",")
}

hgPos <- suppressMessages(
  AnnotationDbi::select(
    x = panOrgDb, keys = hgs,
    columns = c(
      "GID", "mRNA_id",
      "genomeId", "chr", "chr_id", "chr_name", "start", "end", "strand"
    )
  ) 
) %>% 
  dplyr::distinct() %>% 
  dplyr::mutate(
    dplyr::across(.cols = c(start, end), .fns = as.numeric),
    length = end - start + 1
  ) %>% 
  dplyr::filter(genomeId %in% !!genomes) %>% 
  dplyr::group_by(genomeId, chr_id) %>% 
  dplyr::mutate(n = n(), nUniq = length(unique(GID))) %>% 
  # dplyr::add_count(genomeId, chr_id) %>%
  dplyr::filter(n >= 2, nUniq == 2) %>% 
  dplyr::ungroup()

hgGrl <- GenomicRanges::makeGRangesListFromDataFrame(
  df = hgPos, split.field = "genomeId", keep.extra.columns = TRUE
)

hgRange <- unlist(range(hgGrl, ignore.strand = TRUE))
hgRegionGr <- hgRange

# use the inner region boundary instead of outer region boundary
if(opts$inner_region){
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

outCols <- c("regionId", "genomeId", "sampleId", "region", "strand",
             "chr_name", "start", "end")


# optionally, extract HGs and save the column
if(opts$haplotypes){
  hgRegions <- dplyr::rowwise(hgRegions) %>%
    dplyr::mutate(
      hgs = list(
        region_homology_groups(
          pandb = panOrgDb, genome = genomeId, chr = chr_name, start = start,
          end = end, strand = strand, overlapping = opts$overlapping
        )
      )
    ) %>% 
    dplyr::mutate(
      nHgs = length(hgs),
      hgs = paste(hgs, collapse = ";"),
    )
  
  hgRegions %<>% 
    dplyr::ungroup() %>% 
    dplyr::add_count(hgs, name = "grp_n") %>% 
    dplyr::arrange(desc(grp_n)) %>% 
    dplyr::mutate(
      haplotype = paste("haplotype_", dplyr::cur_group_id(), sep = ""),
      .by = hgs
    )
    
  
  outCols <- c(outCols, "hgs", "haplotype")
}

if(is.null(opts$out)){
  dplyr::select(hgRegions, dplyr::all_of(outCols)) %>%
    write.table(sep = "\t", quote = FALSE)
} else {
  dplyr::select(hgRegions, dplyr::all_of(outCols)) %>%
    readr::write_tsv(
      file = paste(opts$dir, "/", opts$out, sep = "")
    )
  
}

