suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(topGO))
suppressPackageStartupMessages(library(AnnotationForge))

## compare homology groups counts between previous and current Pecto pangenomes
rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- panConf$dir
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus) %>%
  dplyr::select(Genome, SpeciesName)

allHgs <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups,
    col_select = 1, col_types = "c"
  )
) %>%
  dplyr::rename(hg_id = "Homology group id") %>% 
  dplyr::mutate(hg_id = paste("hg_", hg_id, sep = ""))

# pangenome GO infor
panDf <- suppressMessages(
  readr::read_tsv(file = panConf$files$go_data, na = "None")
) %>%
  dplyr::mutate(
    GID = paste("hg_", hg_id, sep = "")
  ) %>% 
  dplyr::select(-hg_id)

# chr info
chrInfo <- suppressMessages(
  readr::read_tsv(file = panConf$files$chr_info)
) %>% 
  dplyr::mutate(
    genomeId = paste("g_", genome, sep = "")
  )

# pangenome gene - mRNA - COG - Pfam info
geneInfo <- suppressMessages(readr::read_tsv(panConf$files$gene_info))

geneCoordinates <- dplyr::select(
  geneInfo, mRNA_id, genome, chr_num, start, end, strand, gene_name
  ) %>% 
  dplyr::distinct() %>% 
  dplyr::group_by(genome, chr_num) %>% 
  dplyr::arrange(start, .by_group = TRUE) %>% 
  dplyr::mutate(genePos = 1:n()) %>% 
  dplyr::ungroup()

################################################################################
## combine data to make annotation tables

## homology group - mRNA links
hgDf <- dplyr::select(panDf, GID, mRNA_id, genome, chr) %>%
  dplyr::distinct() %>%
  tidyr::unite(genome, chr, mRNA_id, col = "mRNA_key", sep = "_", remove = FALSE) %>%
  dplyr::select(GID, everything())

hgDf <- dplyr::left_join(
  x = hgDf, y = chrInfo,
  by = c("genome", "chr" = "chr_num")
) %>%
  dplyr::left_join(
    y = geneCoordinates,
    by = c("genome", "chr" = "chr_num", "mRNA_id")
  ) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(genome, chr, start, end), .fns = as.integer
    )
  ) %>% 
  dplyr::distinct()

if (!setequal(hgDf$GID, allHgs$hg_id)) {
  stop("Some missing homology groups")
}

## homology group - GO links
hgGoDf <- dplyr::select(panDf, GID, GO = go_id) %>%
  dplyr::filter(!is.na(GO)) %>%
  dplyr::distinct() %>%
  dplyr::mutate(EVIDENCE = "IEA")

## homology group - COG links
hgCog <- dplyr::select(hgDf, GID, mRNA_id, genome, chr) %>%
  dplyr::left_join(
    y = dplyr::select(geneInfo, mRNA_id, genome, chr_num, starts_with("COG")),
    by = c("genome", "chr" = "chr_num", "mRNA_id")
  ) %>%
  dplyr::select(GID, starts_with("COG")) %>%
  dplyr::distinct()

# homology group - Pfam links
hgPfam <- dplyr::select(hgDf, GID, mRNA_id, genome, chr) %>%
  dplyr::left_join(
    y = dplyr::select(geneInfo, mRNA_id, genome, chr_num, starts_with("pfam_")),
    by = c("genome", "chr" = "chr_num", "mRNA_id")
  ) %>%
  dplyr::select(GID, starts_with("pfam_")) %>%
  dplyr::filter(pfam_id != "None") %>% 
  dplyr::distinct()

## homology group category for pangenome: core/accessory/unique
hgMeta <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups
  )
) %>%
  dplyr::select(GID = "Homology group id", class) %>%
  dplyr::mutate(GID = paste("hg_", GID, sep = ""))

# sed -n '1,/##FASTA/p' interproscan/GCF_024506455.1_ASM2450645v1.interProScan.gff3 | grep -v '^#' | cut -f1 | sort | uniq | wc -l
# sed -n '1,/##FASTA/p' interproscan/GCF_024506455.1_ASM2450645v1.interProScan.gff3 | grep -v '^#' | grep 'GO:' | cut -f1 | sort | uniq | wc -l
################################################################################

## prepare pangenome.db
AnnotationForge::makeOrgPackage(
  hgDf = hgDf,
  hgMeta = hgMeta,
  go = hgGoDf,
  hgCog = hgCog,
  hgPfam = hgPfam,
  version = "1.0.0",
  maintainer = "Lakhansing Pardeshi <lakhansing.pardeshi@wur.nl>",
  author = "Lakhansing Pardeshi",
  outputDir = outDir,
  tax_id = "122277",
  genus = "P",
  species = "ectobacterium.spp.pan",
  goTable = "go",
  verbose = TRUE
)

remove.packages("org.Pectobacterium.spp.pan.eg.db")

devtools::install(
  pkg = file.path(outDir, "org.Pectobacterium.spp.pan.eg.db"),
  upgrade = "never"
)

panOrgDb <- org.Pectobacterium.spp.pan.eg.db
columns(panOrgDb)
keys(panOrgDb)


