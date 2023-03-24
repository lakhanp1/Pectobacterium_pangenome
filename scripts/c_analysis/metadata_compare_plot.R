#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(ggpattern))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(skimr))

## pangenome data summary and comparison

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- confs$analysis$summary$dir

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
pangenomeConf <- confs$data$pangenomes[[pangenome]]

################################################################################

metadata <- suppressMessages(readr::read_tsv(confs$data$reference_data$files$metadata)) %>% 
  dplyr::select(sampleId, source, SpeciesName, geo_loc_country, length, GC_per,
                N50, N90, L50, L90, collection_year)


panMeta <- suppressMessages(readr::read_csv(pangenomeConf$files$metadata)) %>% 
  dplyr::select(Genome, sampleId, virulence, virulence_pcr) %>% 
  dplyr::left_join(y = metadata, by = "sampleId")

metDf <- suppressMessages(readr::read_csv(pangenomeConf$db$metrics$files$per_genome, col_names = T)) %>% 
  dplyr::select(-c(N50, N90, L50, L90)) %>% 
  dplyr::rename_all(
    .funs = ~ stringr::str_replace_all(., "\\s+", "_")
  ) %>% 
  dplyr::left_join(y = panMeta, by = "Genome") %>% 
  dplyr::filter(SpeciesName == "Pectobacterium brasiliense")



metDf %>% 
  # dplyr::filter(! (source %in% c("NCBI"))) %>%
  # dplyr::filter(geo_loc_country != "Netherlands") %>%
  dplyr::filter(virulence == "avirulent") %>%
  dplyr::select(SpeciesName, mRNA_count, length) %>% 
  dplyr::group_by(SpeciesName) %>% 
  skimr::skim()

################################################################################
plotDf <- dplyr::select(
  metDf, Genome, source, sampleId, SpeciesName, geo_loc_country, virulence,
  virulence_pcr, length, GC_per, Median_gene_length, N50, N90, L50, L90,
  Gene_count, mRNA_count, CDS_count, tRNA_count, rRNA_count,
  Homology_groups, Singletons
) %>% 
  dplyr::mutate(
    geo_loc_country = if_else(
      geo_loc_country == "Netherlands", geo_loc_country, "other", missing = "other"
    )
  ) %>% 
  # tidyr::replace_na(replace = list(virulence = "NA", virulence_pcr = "NA")) %>% 
  tidyr::pivot_longer(
    cols = -c(Genome, sampleId, source, SpeciesName, geo_loc_country,
              virulence, virulence_pcr),
    names_to = "field",
    values_to = "value"
  ) %>% 
  dplyr::mutate(
    field = forcats::fct_relevel(
      field, "length", "Gene_count", "mRNA_count", "CDS_count", "tRNA_count"
    )
  )

## geographical location vs mRNA count
pt_loc <- plotDf %>% 
  # dplyr::filter(field %in% c("length")) %>%
  # dplyr::filter(field %in% c("mRNA_count", "length", "N50", "N90", "L50", "L90")) %>%
  dplyr::filter(field %in% c("mRNA_count", "N50")) %>%
  # dplyr::filter(field %in% c("tRNA_count", "rRNA_count")) %>% 
  ggplot2::ggplot(
    mapping = aes(x = geo_loc_country, y = value)
  ) +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0, linewidth = 1) +
  ggbeeswarm::geom_quasirandom(
    mapping = aes(shape = source, color = virulence),
    size = 3, alpha = 0.8
  ) +
  ggpubr::stat_compare_means(method = "t.test") +
  facet_wrap(facets = ~field, nrow = 1, scales = "free_y") +
  labs(
    title = "Pbr isolates from the Netherlands have more protein coding genes"
  ) +
  scale_color_manual(
    values = c("virulent" = "red", "avirulent" = "green"), na.value = "black"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

ggsave(filename = file.path(outDir, "location_vs_ngenes.pdf"), plot = pt_loc, width = 8, height = 8)

## virulence vs mRNA count
pt_vir <- dplyr::filter(plotDf, virulence %in% c("virulent", "avirulent")) %>% 
  dplyr::filter(field %in% c("mRNA_count", "length", "N50")) %>%
  ggplot(mapping = aes(x = virulence, y = value))  +
  geom_boxplot(width = .5, outlier.shape = NA, alpha = 0, linewidth = 1) +
  ggbeeswarm::geom_quasirandom(
    mapping = aes(color = source), size = 3
  ) +
  ggpubr::stat_compare_means(method = "t.test") +
  scale_color_manual(
    values = c("NCBI" = "blue", "NAK" = "green", "NVWA" = "green")
  ) +
  facet_wrap(facets = ~field, nrow = 1, scales = "free_y") +
  labs(
    title = "Virulent Pbr isolates have more protein coding genes than avirulent"
  ) +
  theme_bw(base_size = 14) +
  theme(
    axis.title = element_blank(),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

ggsave(filename = file.path(outDir, "vir_vs_ngenes.pdf"), plot = pt_vir, width = 8, height = 8)
################################################################################

pbrMeta <- dplyr::filter(panMeta, !is.na(collection_year)) %>% 
  dplyr::select(Genome, sampleId, SpeciesName, collection_year, geo_loc_country,
                virulence, virulence_pcr) %>% 
  dplyr::mutate(
    pbr = dplyr::if_else(
      SpeciesName == "Pectobacterium brasiliense", "P. brasiliense", "other"
    ),
    geo_loc_country = dplyr::if_else(
      geo_loc_country == "Netherlands", geo_loc_country, "other", missing = "other"
    ),
    virType = dplyr::case_when(
      virulence == "virulent" & virulence_pcr == "positive" ~ "TP",
      virulence == "virulent" & virulence_pcr == "negative" ~ "FN",
      virulence == "avirulent" & virulence_pcr == "positive" ~ "FP",
      virulence == "avirulent" & virulence_pcr == "negative" ~ "TN",
      TRUE ~ "missing-data"
    ),
    pbr = forcats::fct_relevel(pbr, "other", "P. brasiliense"),
    geo_loc_country = forcats::fct_relevel(geo_loc_country, "Netherlands"),
    virType = forcats::fct_relevel(virType, "missing-data", "TP", "FN", "TN", "FP")
  ) %>% 
  dplyr::filter(SpeciesName == "Pectobacterium brasiliense") %>% 
  dplyr::mutate(
    
  )


pt_pbrTimeline <-  ggplot2::ggplot(pbrMeta) +
  geom_histogram(
    mapping = aes(x = collection_year, fill = forcats::fct_rev(geo_loc_country)),
    binwidth = 2, color = "black"
  ) +
  geom_vline(xintercept = 2015, color = "blue", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Pbr collection after 2015 in the Netherlands"
  ) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(add = c(0, 5))) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = c(0.1, 0.9),
    legend.justification = c(0, 1),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

ggsave(filename = file.path(outDir, "pbr_timeline.pdf"), plot = pt_pbrTimeline, width = 8, height = 6)


pt_fnPbr <- dplyr::filter(pbrMeta, geo_loc_country == "Netherlands") %>% 
  ggplot2::ggplot() +
  ggpattern::geom_histogram_pattern(
    mapping = aes(x = collection_year, fill = virType, pattern_shape = virType, pattern_fill = virType),
    binwidth = 1, pattern = "regular_polygon", pattern_density = 0.6, pattern_grid = "hex",
    color = "black", pattern_spacing = 0.02, pattern_colour = "black"
  ) +
  scale_fill_manual(
    values = c("TP" = "red", "FP" = alpha("green", 0.6),
               "TN" = "green", "FN" = alpha("red", 0.6), "missing-data" = "grey")
  ) +
  scale_pattern_shape_manual(
    values = c("TP" = "null", "FP" = "convex6", "TN" = "null", "FN" = "convex6", "missing-data" = "null")
  ) +
  scale_pattern_fill_manual(
    values = c("TP" = "red", "FP" = "red", "TN" = "green", "FN" = "green", "missing-data" = "grey")
  ) +
  labs(
    title = "FN-Pbr emergence in the Netherlands"
  ) +
  scale_x_continuous(expand = expansion(add = 0)) +
  scale_y_continuous(expand = expansion(add = c(0, 5))) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = c(0.1, 0.9),
    legend.justification = c(0, 1),
    legend.key.size = unit(1, "cm"),
    legend.text = element_text(size = 20),
    legend.title = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )

pt_fnPbr

ggsave(filename = file.path(outDir, "fn_pbr_NL.pdf"), plot = pt_fnPbr, width = 8, height = 6)



