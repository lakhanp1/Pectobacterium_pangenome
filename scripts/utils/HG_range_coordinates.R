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
  help = "[optional] a COMMA separated list of genome ids. Default: use all genomes"
)

parser <- optparse::add_option(
  parser, default = FALSE,
  opt_str = c("--inner_region"), type = "logical", action = "store_true",
  help = "LOGICAL: TRUE: use the inner boundry of the region, othewise outer"
)

parser <- optparse::add_option(
  parser, default = NULL,
  opt_str = c("--max_genes"), type = "integer", action = "store",
  help = "INTEGER: Maximum number of genes expected for a region"
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
  parser, default = NULL,
  opt_str = c("-o", "--out"), type = "character", action = "store",
  help = "Output file name"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$hgs))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# opts$hgs <- "hg_22427625,hg_22427622"
# opts$out <- "analysis/pangenome_v2/carotovoricin/ctv_nitric_oxide/hg_regions.tab"
# opts$inner_region <- FALSE
# opts$haplotypes <- TRUE
# opts$max_genes <- 10
# opts$genomes <- "g_162,g_85,g_418,g_294,g_20,g_116,g_64,g_1,g_10,g_100,g_262,g_149"
################################################################################

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

hgs <- stringr::str_split_1(string = opts$hgs, pattern = ",")

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

if(!dir.exists(dirname(opts$out))){
  stop(
    paste("Output directory does not exist:", dirname(opts$out)),
    call. = FALSE
  )
}

if(is.null(opts$max_genes)){
  opts$max_genes <- Inf
}

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

# ensure that there is only the closest hg1-hg2 pair in the output
hg1Gr <- GenomicRanges::makeGRangesFromDataFrame(
  df = dplyr::filter(hgPos, GID == !!hgs[1]),
  seqnames.field = "chr_id",
  keep.extra.columns = TRUE
)

hg2Gr <- GenomicRanges::makeGRangesFromDataFrame(
  df = dplyr::filter(hgPos, GID == !!hgs[2]),
  seqnames.field = "chr_id",
  keep.extra.columns = TRUE
)

hg1ToHg2 <- GenomicRanges::nearest(
  x = hg1Gr, subject = hg2Gr, ignore.strand = TRUE
)

hg2ToHg1 <- GenomicRanges::nearest(
  x = hg2Gr, subject = hg1Gr, ignore.strand = TRUE
)

hgGrPairs <- tibble::tibble(
  hg1 = 1:length(hg1ToHg2),
  hg2 = hg1ToHg2
) %>%
  dplyr::bind_rows(tibble::tibble(hg1 = hg2ToHg1, hg2 = 1:length(hg2ToHg1))) %>%
  dplyr::distinct()

hgGrPairs$dist <- GenomicRanges::distance(
  x = hg1Gr[hgGrPairs$hg1], y = hg2Gr[hgGrPairs$hg2], ignore.strand = TRUE
)

hgGrPairs$chr <- as.character(GenomicRanges::seqnames(hg1Gr[hgGrPairs$hg1]))

closestGrPairs <- dplyr::group_by(hgGrPairs, chr) %>%
  dplyr::arrange(dist, .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

# new closest hg1-hg2 pairs
hg1GrNew <- hg1Gr[closestGrPairs$hg1]
hg2GrNew <- hg2Gr[closestGrPairs$hg2]

# convert to GRangesList for region border extraction
hgGrl <- GenomicRanges::makeGRangesListFromDataFrame(
  df = as.data.frame(c(hg1GrNew, hg2GrNew)), split.field = "genomeId",
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames"
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
  dplyr::select(genomeId, chr_id = seqnames, start, end) %>%
  tibble::as_tibble() %>%
  dplyr::left_join(
    y = as.data.frame(hg1GrNew) %>%
      dplyr::select(genomeId, chr_id = seqnames, chr, chr_name, strand) %>%
      dplyr::distinct(),
    by = c("genomeId", "chr_id")
  ) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, genomeId, sampleId, SpeciesName),
    by = "genomeId"
  ) %>%
  dplyr::mutate(
    region = paste(genomeId, "#", chr_name, ":", start, "-", end, sep = ""),
    SpeciesName = stringr::str_replace_all(
      string = SpeciesName, pattern = "\\.?\\s+", replacement = "_"
    ),
    regionId = paste(genomeId, ".", SpeciesName, sep = ""),
    length = abs(end - start)
  ) %>%
  dplyr::select(
    regionId, genomeId, sampleId, chr_name, start, end, strand, region, length
  )

outCols <- c("regionId", "genomeId", "sampleId", "region", "strand",
             "chr_name", "start", "end", "length")


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

  hgRegions <- dplyr::ungroup(hgRegions) %>%
    dplyr::add_count(hgs, name = "grp_n") %>%
    dplyr::arrange(desc(grp_n)) %>%
    dplyr::mutate(
      haplotype = paste("haplotype_", dplyr::cur_group_id(), sep = ""),
      .by = hgs
    ) %>%
    dplyr::filter(nHgs <= opts$max_genes)

  outCols <- c(outCols, "nHgs", "hgs", "haplotype")
}



if(is.null(opts$out)){
  dplyr::select(hgRegions, dplyr::all_of(outCols)) %>%
    write.table(sep = "\t", quote = FALSE, row.names = FALSE)
} else {
  dplyr::select(hgRegions, dplyr::all_of(outCols)) %>%
    readr::write_tsv(file = opts$out)

}

