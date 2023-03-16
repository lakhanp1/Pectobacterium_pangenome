#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(skimr))

## compare homology groups counts between previous and current Pecto pangenomes
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

grpStats <- suppressMessages(readr::read_tsv(file.path(outDir, "group_stats.tab"))) %>% 
  dplyr::mutate(
    data = forcats::as_factor(data),
    group = forcats::as_factor(group)
  )

pt_hgStats <- ggplot2::ggplot(grpStats) +
  geom_bar(
    mapping = aes(x = data, y = count, fill = forcats::fct_rev(group)),
    stat = "identity",
    position = position_stack()
  ) +
  scale_y_continuous(expand = expansion(add = c(0, 100))) +
  scale_fill_manual(
    name = NULL,
    values = c("core" = confs$colors$core, "accessory" = confs$colors$accessory,
               "unique" = confs$colors$unique)
    ) +
  theme_bw(base_size = 24) +
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title = element_blank(),
    axis.text.x = element_text(face = "italic")
  )

ggsave(filename = file.path(outDir, "hg_cmp_stats.pdf"), plot = pt_hgStats, width = 6, height = 6)

