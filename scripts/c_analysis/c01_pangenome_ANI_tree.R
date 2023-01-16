suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))

## generate UPGMA and NJ tree using ANI for pangenome

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

outDir <- confs$analysis$phylogeny$dir
genus <- "Pectobacterium"

pt_theme <- theme_bw(base_size = 14) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, face = "bold")
  )

################################################################################
genusPattern <- paste("(", genus, " )(?!sp\\.)", sep = "")

sampleInfo <- suppressMessages(
  readr::read_csv(file = confs$data$pangenomes$pectobacterium.v2$files$metadata)
) %>% 
  dplyr::mutate(
    SpeciesName = stringi::stri_replace(
      str = SpeciesName, regex = genusPattern, replacement = "P. "
    ),
    SpeciesName = stringi::stri_replace(
      str = SpeciesName, regex = "((\\w)[^ ]+ )((\\w)[^ ]+ )(subsp\\..*)",
      replacement = "$2. $4. $5"
    ),
    nodeLabs = stringr::str_c(sampleName, " (", SpeciesName,")", sep = "")
  )

sampleInfoList <- dplyr::select(
  sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
) %>% 
  purrr::transpose() %>% 
  purrr::set_names(nm = purrr::map(., "sampleId"))

genomeIds <- dplyr::pull(sampleInfo, Genome, name = sampleId)

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


distMat <- tidyr::pivot_wider(
  data = aniDf,
  id_cols = "g1",
  names_from = "g2",
  values_from = "dist"
) %>% 
  tibble::column_to_rownames(var = "g1") %>% 
  as.matrix() %>% 
  as.dist()

aniMat <- tidyr::pivot_wider(
  data = aniDf,
  id_cols = "g1",
  names_from = "g2",
  values_from = "ani"
) %>% 
  tibble::column_to_rownames(var = "g1") %>% 
  as.matrix()

if(!all(rownames(as.matrix(distMat)) == rownames(aniMat))){
  stop("rownames did not match")
}
################################################################################

# plot(hclust(distMat))
treeUpgma <- ape::as.phylo(hclust(d = distMat, method = "average")) %>% 
  ape::ladderize()

ape::write.tree(
  phy = treeUpgma, tree.names = "ANI_UPGMA",
  file = confs$analysis$phylogeny$files$tree_ani_upgma
)

treeNj <- ape::nj(distMat) %>% 
  ape::ladderize()

ape::write.tree(
  phy = treeNj, tree.names = "ANI_NJ",
  file = confs$analysis$phylogeny$files$tree_ani_nj
)






