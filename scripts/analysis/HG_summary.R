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
suppressPackageStartupMessages(library(matrixStats))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

## 1) process homology groups to transpose the table where
## rows are genomes and columns are homology group identifiers
## 2) plot heatmap of homology group across pangenome

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/GO_enrichment/enrichment_functions.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

analysisName <- "homology_groups"
outDir <- confs$analysis$homology_groups$dir
outPrefix <- file.path(outDir, analysisName)
panOrgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

spOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$core_snp_ml$files$species_order)
)

hgTable <- AnnotationDbi::select(
  x = panOrgDb, keys = keys(panOrgDb),
  columns = c("GID", "genomeId", "class")
) %>%
  dplyr::rename(hgId = GID) %>%
  dplyr::count(hgId, genomeId, name = "nGenes")

hgMeta <- AnnotationDbi::select(
  x = panOrgDb, keys = keys(panOrgDb),
  columns = c("GID", "class")
) %>%
  dplyr::rename(hgId = GID) %>%
  dplyr::mutate(
    class = dplyr::if_else(
      condition = class == "core & single copy orthologous",
      true = "core", false = class
    )
  )

################################################################################
## GO enrichment for homology groups

spNames <- dplyr::count(sampleInfo, SpeciesName) %>%
  # dplyr::filter(n >= 20) %>%
  dplyr::pull(SpeciesName)

sppGrpStats <- NULL
sppGrpGo <- NULL

## get species wise core, accessory, unique group stats and GO
for (sp in spNames) {
  spGenomes <- dplyr::filter(sampleInfo, SpeciesName == .env$sp) %>%
    dplyr::pull(genomeId)
  
  cat(sp, ": ", length(spGenomes), "\n")
  
  hgSum <- sub_pangenome_hgs(pandb = panOrgDb, genomes = spGenomes) %>% 
    dplyr::left_join(
      y = dplyr::rename(hgMeta, pangenome_class = class),
      by = "hgId"
    )
  
  ## group stats
  sppGrpStats <- dplyr::count(hgSum, pangenome_class, class, name = "count") %>%
    dplyr::mutate(SpeciesName = .env$sp, nHgs = !!nrow(hgSum)) %>%
    dplyr::bind_rows(sppGrpStats)
  
  ## GO enrichment
  sppGrpGo <- dplyr::group_by(hgSum, class) %>%
    dplyr::group_modify(
      .f = ~ topGO_enrichment(genes = .x$hgId, orgdb = panOrgDb)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(SpeciesName = .env$sp) %>%
    dplyr::bind_rows(sppGrpGo)
}

################################################################################
## pangenome wide GO enrichment
panGo <- hgMeta  %>%
  dplyr::group_by(class) %>%
  dplyr::group_modify(
    .f = ~ topGO_enrichment(genes = .x$hgId, orgdb = panOrgDb)
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SpeciesName = "pangenome")

sppGrpGo <- dplyr::bind_rows(panGo, sppGrpGo) %>%
  dplyr::select(SpeciesName, everything())

readr::write_tsv(
  x = sppGrpGo, file = confs$analysis$homology_groups$files$spp_group_go
)

################################################################################

ptDf <- dplyr::group_by(panGo, class) %>%
  dplyr::arrange(pvalue, .by_group = TRUE) %>%
  dplyr::slice(1:8) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    class = forcats::fct_relevel(class, "core", "accessory", "unique")
  )

pt_go <- enrichment_bar(df = ptDf, title = "pangenome GO", colorCol = "class")

(pt_go <- pt_go +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  scale_x_discrete(labels = label_wrap(60)) +
  scale_fill_manual(
    values = c(
      "core" = confs$colors$core, "accessory" = confs$colors$accessory,
      "unique" = confs$colors$unique
    )
  )
)

ggsave(
  plot = pt_go, filename = file.path(outDir, "pangenome_GO.pdf"),
  width = 8, height = 6
)


################################################################################
# core-accessory-unique stats for each species
hgStats <- dplyr::count(hgMeta, class, name = "count") %>%
  dplyr::mutate(
    SpeciesName = "Pangenome",
    nHgs = sum(count),
    pangenome_class = class
  ) %>%
  dplyr::bind_rows(sppGrpStats) %>%
  dplyr::select(SpeciesName, pangenome_class, class, count, nHgs) %>% 
  dplyr::mutate(
    fraction = round(count / nHgs, digits = 3),
    class = forcats::fct_relevel(class, "core", "accessory", "unique"),
    SpeciesName = forcats::fct_relevel(SpeciesName, "Pangenome", !!!spOrder$SpeciesName)
  )

readr::write_tsv(
  x = hgStats, file = confs$analysis$homology_groups$files$spp_group_stats
)

species_to_show <- dplyr::count(sampleInfo, SpeciesName) %>%
  dplyr::filter(n >= 20) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::mutate(
    SpeciesName = forcats::as_factor(SpeciesName)
  )

pt_stats <- dplyr::left_join(
  x = species_to_show, y = hgStats, by = "SpeciesName"
) %>% 
  ggplot() +
  geom_bar(
    mapping = aes(
      x = count, y = SpeciesName,
      fill = forcats::fct_rev(class)
    ),
    stat = "sum", position = position_dodge(), width = 0.8
  ) +
  scale_fill_manual(
    values = c(
      confs$colors$core, confs$colors$accessory, confs$colors$unique
    ),
    breaks = c("core", "accessory", "unique")
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold.italic"),
    legend.position = "bottom",
    legend.title = element_blank()
  )


ggsave(
  plot = pt_stats, filename = file.path(outDir, "pangenome_spp_stats.pdf"),
  width = 6, height = 4
)

################################################################################
## plot Heap's law alpha for multiple species in pangenome
heaps <- suppressMessages(readr::read_tsv(file.path(outDir, "heaps_law.tab"))) %>%
  dplyr::filter(species != "pangenome") %>%
  dplyr::mutate(
    complete = "complete",
    species = forcats::fct_relevel(species, !!!levels(species_to_show$SpeciesName))
  )

(pt_alpha <- ggplot(data = heaps) +
    geom_bar(
      mapping = aes(y = species, x = alpha),
      fill = "black", stat = "identity", width = 0.8
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.01))) +
    theme_bw(base_size = 16) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_text(face = "bold"),
      axis.text.y = element_text(face = "bold.italic")
    )
)

ggsave(
  plot = pt_alpha, filename = file.path(outDir, "pangenome_spp_alpha.pdf"),
  width = 5, height = 4
)

################################################################################
