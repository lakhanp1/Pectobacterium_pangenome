suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))

## intra-species and inter-species ANI score distribution

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
outGroup <- confs$analysis$phylogeny$outgroup
outDir <- confs$analysis$ANI$dir


################################################################################
sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

genomeIds <- dplyr::pull(sampleInfo, Genome, name = sampleId)

speciesOrd <- names(sort(table(sampleInfo$SpeciesName), decreasing = T))

aniDf <- suppressMessages(readr::read_tsv(
  file = confs$analysis$ANI$files$fastani_out,
  col_names = c("id1", "id2", "ani", "mapped", "total")
))

aniDf %<>% dplyr::mutate(
  dplyr::across(
    .cols = c(id1, id2),
    .fns = ~stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1"),
  ),
  dist = 1 - (ani/100)
) %>% 
  dplyr::mutate(
    g1 = genomeIds[id1],
    g2 = genomeIds[id2]
  ) %>% 
  dplyr::filter(!is.na(g1) | !is.na(g2)) %>% 
  dplyr::arrange(g1, g2)


aniMat <- tidyr::pivot_wider(
  data = aniDf,
  id_cols = "g1",
  names_from = "g2",
  values_from = "ani"
) %>% 
  tibble::column_to_rownames(var = "g1") %>% 
  as.matrix()

aniMat[upper.tri(aniMat)] <- NA

intraSpAni <- as.data.frame(aniMat) %>% 
  tibble::rownames_to_column(var = "g1") %>% 
  tidyr::pivot_longer(
    cols = -g1,
    names_to = "g2",
    values_to = "ani"
  ) %>% 
  dplyr::filter(!is.na(ani), g1 != g2) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, sp1 = SpeciesName, type1 = type_material),
    by = c("g1" = "Genome")
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, sp2 = SpeciesName, type2 = type_material),
    by = c("g2" = "Genome")
  ) %>% 
  dplyr::filter(sp1 == sp2) %>% 
  dplyr::mutate(
    aniType = "intra-species",
    species = sp1
  )

typeStrainAni <- dplyr::select(aniDf, g1, g2, ani) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, sp1 = SpeciesName, type1 = type_material),
    by = c("g1" = "Genome")
  ) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, sp2 = SpeciesName, type2 = type_material),
    by = c("g2" = "Genome")
  ) %>% 
  dplyr::filter(g1 != g2, sp1 != sp2, !is.na(type2)) %>% 
  dplyr::mutate(
    aniType = "inter-species",
    species = sp1
  )

aniCmp <- dplyr::bind_rows(intraSpAni, typeStrainAni) %>% 
  dplyr::mutate(
    species = forcats::fct_relevel(species, !!!speciesOrd),
    aniType = forcats::fct_relevel(aniType, "intra-species", "inter-species")
  )

################################################################################

pt_aniDistr <- ggplot2::ggplot(
  data = aniCmp,
  mapping = aes(x = aniType, y = ani, color = aniType)
) +
  geom_hline(yintercept = 95, color = "black", linetype = "dashed") +
  ggbeeswarm::geom_quasirandom() +
  facet_wrap(facets = vars(species), scales = "free_y", nrow = 4) +
  scale_color_manual(
    name = NULL,
    values = c("intra-species" = "#FE5F55", "inter-species" = "#227C9D")
  ) +
  scale_y_continuous(expand = expansion(add = 1), limits = c(min(aniCmp$ani), 100)) +
  labs(y = "ANI") +
  guides(colour = guide_legend(override.aes = list(size=5))) +
  theme_bw(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    legend.position = "bottom",
    axis.title.x = element_blank(),
    axis.text.x.bottom = element_blank(),
    legend.text = element_text(size = 16), 
    strip.text.x = element_text(face = "italic")
  )

pt_aniDistr

ggsave(
  filename = file.path(outDir, "species_ANI_distribution.png"),
  plot = pt_aniDistr, width = 16, height = 12
)




