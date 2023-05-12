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
orgDb <- org.Pectobacterium.spp.pan.eg.db
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

spOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny$ani_upgma$files$species_order)
)

hgs <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups
  )
) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "Genome_", "")) %>% 
  dplyr::rename_with(.fn = ~tolower(.x)) %>% 
  dplyr::mutate(homology_group_id = as.character(homology_group_id))


hgTable <- dplyr::select(hgs, hg = homology_group_id, !!!sampleInfo$Genome) %>% 
  tidyr::pivot_longer(
    cols = -hg, names_to = "Genome", values_to = "nGenes"
  )

g2hg <- tidyr::pivot_wider(
  hgTable,
  id_cols = Genome, names_from = hg, values_from = nGenes
)

readr::write_tsv(
  x = g2hg, file = confs$analysis$homology_groups$files$groups
)

hgMeta <- dplyr::select(hgs, homology_group_id, class)

readr::write_tsv(x = hgMeta, file = confs$analysis$homology_groups$files$groups_meta)
################################################################################
## GO enrichment for homology groups

## binary matrix for homology_group PAV
hgBinaryMat <- dplyr::mutate(
  hgTable,
  nGenes = dplyr::if_else(nGenes == 0, 0, 1)
) %>% 
  tidyr::pivot_wider(
    id_cols = Genome, names_from = hg, values_from = nGenes
  ) %>% 
  tibble::column_to_rownames(var = "Genome") %>% 
  as.matrix()

spNames <- dplyr::count(sampleInfo, SpeciesName) %>% 
  dplyr::filter(n >= 20) %>% 
  dplyr::pull(SpeciesName)

sppGrpStats <- NULL
sppGrpGo <- NULL

## get species wise core, accessory, unique group stats and GO
for (sp in spNames) {
  
  spGenomes <- dplyr::filter(sampleInfo, SpeciesName == .env$sp) %>% 
    dplyr::pull(Genome)
  
  cat(sp, length(spGenomes), "\n")
  
  hgSum <- matrixStats::colSums2(
    x = hgBinaryMat, useNames = T,
    rows = which(rownames(hgBinaryMat) %in% spGenomes)
  ) %>% 
    tibble::enframe(name = "hg", value = "nGenomes") %>% 
    dplyr::filter(nGenomes != 0) %>% 
    dplyr::mutate(
      class = dplyr::case_when(
        nGenomes == 1 ~ "unique",
        nGenomes == !!length(spGenomes) ~ "core",
        nGenomes < !!length(spGenomes) & nGenomes > 1 ~ "accessory"
      )
    )
  
  ## group stats
  sppGrpStats <- dplyr::count(hgSum, class, name = "count") %>% 
    dplyr::mutate(SpeciesName = .env$sp, fraction = count/!!nrow(hgSum)) %>% 
    dplyr::bind_rows(sppGrpStats)
  
  ## GO enrichment 
  sppGrpGo <- dplyr::group_by(hgSum, class) %>% 
    dplyr::group_modify(
      .f = ~topGO_enrichment(genes = .x$hg, orgdb = orgDb)
    ) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(SpeciesName = .env$sp) %>% 
    dplyr::bind_rows(sppGrpGo)
  
}


################################################################################
## pangenome wide GO enrichment
panGo <- hgMeta %>% 
  dplyr::mutate(
    class = dplyr::if_else(
      condition = class == "core & single copy orthologous",
      true = "core", false = class
    )
  ) %>% 
  dplyr::group_by(class) %>% 
  dplyr::group_modify(
    .f = ~topGO_enrichment(genes = .x$homology_group_id, orgdb = orgDb)
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(SpeciesName = "pangenome")

sppGrpGo <- dplyr::bind_rows(panGo, sppGrpGo) %>% 
  dplyr::select(SpeciesName, everything())

readr::write_tsv(
  x = sppGrpStats, file = confs$analysis$homology_groups$files$spp_group_stats
)

readr::write_tsv(
  x = sppGrpGo, file = confs$analysis$homology_groups$files$spp_group_go
)

################################################################################

ptDf <- dplyr::group_by(panGo, class) %>% 
  dplyr::arrange(pvalue, .by_group = TRUE) %>% 
  dplyr::slice(1:7) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(
    class = forcats::fct_relevel(class, "core", "accessory", "unique")
  )

pt_go <- enrichment_bar(df = ptDf, title = "pangenome GO", colorCol = "class")

(pt_go <- pt_go +
    scale_y_continuous(expand = expansion(mult = c(0,0.2))) +
    scale_x_discrete(labels = label_wrap(60)) +
    scale_fill_manual(
      values = c("core" = confs$colors$core, "accessory" = confs$colors$accessory,
                 "unique" = confs$colors$unique)
    )
)

ggsave(
  plot = pt_go, filename = file.path(outDir, "pangenome_GO.pdf"),
  width = 8, height = 6
)


################################################################################
# core-accessory-unique stats for each species
sppGrpStats <-
  # dplyr::select(panGo, class, count = inputSize) %>% 
  # dplyr::distinct() %>% 
  # dplyr::mutate(
  #   SpeciesName = "Pangenome",
  #   fraction = count / sum(count)
  # ) %>% 
  # dplyr::bind_rows(
  suppressMessages(
    readr::read_tsv(confs$analysis$homology_groups$files$spp_group_stats)
    # )
  ) %>% 
  dplyr::mutate(
    class = forcats::fct_relevel(class, "core", "accessory", "unique"),
    SpeciesName = forcats::fct_relevel(SpeciesName, "Pangenome", !!!spOrder$SpeciesName)
  )

(pt_stats <- ggplot(data = sppGrpStats) +
    geom_bar(
      mapping = aes(
        x = count, y = forcats::fct_rev(SpeciesName),
        fill = forcats::fct_rev(class)
      ),
      stat = "identity", position = position_dodge(), width = 0.8
    ) +
    scale_fill_manual(
      values = c("core" = confs$colors$core, "accessory" = confs$colors$accessory,
                 "unique" = confs$colors$unique)
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
)

ggsave(
  plot = pt_stats, filename = file.path(outDir, "pangenome_spp_stats.pdf"),
  width = 5, height = 4
)

################################################################################
## plot Heap's law alpha for multiple species in pangenome
heaps <- suppressMessages(readr::read_tsv(file.path(outDir, "heaps_law.tab"))) %>% 
  dplyr::filter(species != "pangenome") %>% 
  dplyr::mutate(
    complete = "complete",
    species = forcats::fct_relevel(species, !!!spOrder$SpeciesName)
  )

(pt_alpha <- ggplot(data = heaps) +
    geom_bar(
      mapping = aes(y = forcats::fct_rev(species), x = alpha),
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


