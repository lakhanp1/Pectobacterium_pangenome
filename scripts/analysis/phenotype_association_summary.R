#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(viridisLite))
suppressPackageStartupMessages(library(openxlsx))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

## plot PAV for the homology groups specific for a particular phenotype

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/GO_enrichment/enrichment_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/association_analysis.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

phenotype <- "assay_FN"
treeMethod <- "ani_upgma" # ani_upgma, kmer_nj
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

outDir <- file.path(confs$analysis$association$path, phenotype)
outPrefix <- file.path(outDir, phenotype)

## sequence info file for mRNAs across all genomes belonging to homology groups
## specific to a phenotype of interest
file_associatedSeqInfo <- paste(outPrefix, ".pheno_specific.seq_info.txt", sep = "")
orgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

## genomes that are used for phenotype association analysis
associatedGenomes <- get_phenotype_association_genomes(
  phenotype = phenotype,
  confFile = panConf$analysis_confs$files$clade_association
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>% dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!speciesOrder$SpeciesName)
)

## get phenotype specific groups
specGrps <- phenotype_specific_groups(
  phenotype = phenotype, panConf = panConf
) %>% as.character()


## get the best genome
hgSeqInfo <- AnnotationDbi::select(
  x = orgDb,
  keys = specGrps,
  columns = c(
    "genome", "chr", "chr_id", "chr_name", "start", "end", "strand", "mRNA_id", "mRNA_key"
  )
) %>%
  dplyr::rename(Genome = genome, hg_id = GID) %>%
  dplyr::mutate(
    dplyr::across(.cols = c(chr, start, end), .fns = as.integer)
  ) %>%
  tibble::as_tibble()


genomeSummary <- dplyr::add_count(hgSeqInfo, Genome, name = "nHg") %>%
  dplyr::group_by(Genome, nHg) %>%
  dplyr::summarise(nChr = length(unique(chr)), .groups = "drop") %>%
  dplyr::left_join(
    dplyr::select(
      sampleInfo, Genome, SpeciesName, geo_loc_country, host, isolation_source,
      env_broad_scale, collection_year
    ),
    by = "Genome"
  ) %>%
  dplyr::arrange(desc(nHg), nChr)

readr::write_tsv(
  x = genomeSummary,
  file = paste(outPrefix, ".specific_hgs.pangenome_summary.tab", sep = "")
)

bestGenome <- associatedGenomes$compare[1]

hgAnnotation <- dplyr::filter(hgSeqInfo, Genome == !!bestGenome) %>%
  dplyr::arrange(chr, start) %>%
  dplyr::mutate(
    label = paste(mRNA_id, "| ", chr, ":", start, sep = ""),
    hg_group = chr_id
  )

################################################################################
# functional annotation for homology groups
# COG
grpCog <- AnnotationDbi::select(
  x = orgDb,
  keys = hgAnnotation$hg_id,
  columns = c("GID", "COG_description", "COG_category")
) %>%
  dplyr::rename(hg_id = GID) %>%
  dplyr::left_join(y = hgAnnotation, by = "hg_id") %>%
  dplyr::select(
    hg_id, Genome, chr, start, end, strand, mRNA_id,
    starts_with("COG")
  )

# GO assignment
goMap <- AnnotationDbi::select(
  x = orgDb,
  keys = hgAnnotation$hg_id,
  columns = c("GID", "GO")
)

goDesc <- AnnotationDbi::select(
  x = GO.db, keys = na.omit(goMap$GO), column = "TERM", keytype = "GOID"
) %>%
  dplyr::distinct() %>%
  dplyr::mutate(description = paste("(", GOID, ") ", TERM, sep = "")) %>%
  dplyr::select(-TERM) %>%
  dplyr::left_join(y = goMap, by = c("GOID" = "GO")) %>%
  dplyr::group_by(GID) %>%
  dplyr::summarise(GO = paste(description, collapse = "; "), .groups = "drop")

goMap <- dplyr::left_join(x = hgAnnotation, y = goDesc, by = c("hg_id" = "GID")) %>%
  dplyr::select(
    hg_id, Genome, chr, start, end, strand, mRNA_id, GO
  )

# topGO
topGoDf <- topGO_enrichment(genes = specGrps, orgdb = orgDb)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb = wb, sheetName = "COG")
openxlsx::writeData(
  wb = wb, sheet = 1, x = grpCog,
  startCol = 1, startRow = 1, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 2, firstActiveCol = 2)

openxlsx::addWorksheet(wb = wb, sheetName = "GO")
openxlsx::writeData(
  wb = wb, sheet = 2, x = goMap,
  startCol = 1, startRow = 1, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 2, firstActiveRow = 2, firstActiveCol = 2)

openxlsx::addWorksheet(wb = wb, sheetName = "topGO")
openxlsx::writeData(
  wb = wb, sheet = 3, x = topGoDf,
  startCol = 1, startRow = 1, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 3, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(wb = wb, file = paste(outPrefix, ".functions.xlsx", sep = ""), overwrite = TRUE)

################################################################################
# prepare homology group PAV matrix from pan.db
hgMat <- homology_groups_mat(pandb = orgDb, type = "cnv", groups = hgAnnotation$hg_id)

hgMat <- hgMat[rawTree$tip.label, ]

htList <- homology_group_heatmap(
  mat = hgMat, phy = rawTree, metadata = sampleInfo,
  hgAn = hgAnnotation, width = c(9, 20),
  markGenomes = associatedGenomes
)

pdf(file = paste(outPrefix, ".pheno_hg_association.pdf", sep = ""), width = 14, height = 8)
ComplexHeatmap::draw(
  object = htList,
  main_heatmap = "hg",
  row_dend_side = "left",
)
dev.off()

################################################################################
