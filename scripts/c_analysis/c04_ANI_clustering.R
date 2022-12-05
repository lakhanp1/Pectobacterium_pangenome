suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggpubr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(ggtreeExtra))
suppressPackageStartupMessages(library(ggstar))

## ANI distance matrix clustering using UPGMA to
## 1) verify species clustering
## 2) estimate the outgroup for rooting of phylogenetic tree

rm(list = ls())

################################################################################
set.seed(124)

file_fastani <- here::here("analysis", "02_fastANI", "ANI_results")
file_metadata <- here::here("data/pangenomes/pectobacterium.v2", "genomes_metadata.csv")
file_sampleSummary <- here::here("data/reference_data", "sample_metadata.tsv")

outDir <- here::here("analysis", "02_fastANI")

pt_theme <- theme_bw() +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 14, face = "bold"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 15, face = "bold")
  )

################################################################################
metadata <- suppressMessages(readr::read_tsv(file = file_sampleSummary)) %>% 
  dplyr::select(sampleId, length, GC_per, n_contigs, N50,L50)

sampleInfo <- suppressMessages(readr::read_csv(file = file_metadata)) %>% 
  dplyr::mutate(
    SpeciesName = stringr::str_replace(
      string = SpeciesName, pattern = "Pectobacterium", replacement = "P."
    ),
    nodeLabs = stringr::str_c(sampleName, " (", SpeciesName,")", sep = "")
  ) %>% 
  dplyr::select(id, everything()) %>% 
  dplyr::left_join(y = metadata, by = c("id" = "sampleId"))


sampleInfoList <- dplyr::select(
  sampleInfo, id, sampleName, SpeciesName, strain, nodeLabs
) %>% 
  purrr::transpose() %>% 
  purrr::set_names(nm = purrr::map(., "id"))

aniDf <- suppressMessages(readr::read_tsv(
  file = file_fastani,
  col_names = c("g1", "g2", "ani", "mapped", "total")
))

aniDf %<>% dplyr::mutate(
  dplyr::across(
    .cols = c(g1, g2),
    .fns = ~stringr::str_replace(string = .x, pattern = ".*/(.*).fna", replacement = "\\1"),
  ),
  dist = 1 - (ani/100)
) %>% 
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

################################################################################

# plot(hclust(distMat))
treeUpgma <- ape::as.phylo(hclust(d = distMat, method = "average"))
ape::write.tree(
  phy = treeUpgma, tree.names = "ANI_UPGMA",
  file = file.path(outDir, "ANI_UPGMA.newick")
)

treeNj <- ape::nj(distMat)

corr_nj <- cor(as.vector(distMat), as.vector(as.dist(cophenetic(treeNj))))
corr_upgma <- cor(as.vector(distMat), as.vector(as.dist(cophenetic(treeUpgma))))

## cophenetic distance calculation for NJ and UPGMA tree
copheneticDist <- tibble::tibble(
  dist = as.vector(distMat),
  UPGMA = as.vector(as.dist(cophenetic(treeUpgma))),
  NJ =  as.vector(as.dist(cophenetic(treeNj)))
) %>% 
  tidyr::pivot_longer(cols = -dist, names_to = "tree", values_to = "cophenetic") %>% 
  dplyr::mutate(tree = forcats::as_factor(tree))

pt_copheneticCorr <- ggplot2::ggplot(
  data = copheneticDist,
  mapping = aes(x = dist, y = cophenetic, group = tree)
) +
  geom_point(alpha = 0.7) +
  ggpubr::stat_cor(size = 8) +
  labs(
    title = "Cophenetic distance for UPGMA and NJ trees build using ANI",
    x = "Cophenetic distance", y = "1 - (ANI/100)"
  ) +
  facet_grid(cols = vars(tree)) +
  pt_theme +
  theme(strip.text = element_text(size = 14, face = "bold"))

ggsave(
  plot = pt_copheneticCorr, filename = file.path(outDir, "copheneticDistCorr.png"),
  height = 6, width = 10
)

################################################################################
aniTreeTable <- as_tibble(treeUpgma)

# treeUpgma <- as_tibble(treeUpgma) %>% 
#   dplyr::full_join(y = sampleInfo, by = c("label" = "id")) %>% 
#   treeio::as.treedata()
# 
# treeNj <- as_tibble(treeNj) %>% 
#   dplyr::full_join(y = sampleInfo, by = c("label" = "id")) %>% 
#   treeio::as.treedata()

## Genome, AssemblyAccession, AssemblyName, SpeciesName, , type_material,
## taxonomy_check_status, strain, virulence, virulence_pcr,
## geo_loc_country, host, isolation_source, collected_by, env_broad_scale
pt_upgma1 <- ggtree::ggtree(tr = treeUpgma) +
  labs(title = "UPGMA tree")

pt_upgma2 <- pt_upgma1 %<+% sampleInfo + 
  geom_point(
    mapping = aes(shape = SpeciesName, color = SpeciesName),
    size = 4
  ) +
  scale_shape_manual(name = "outgroup", values = c("P. fontis" = 16)) +
  scale_color_manual(name = "outgroup", values = c("P. fontis" = "red")) +
  ggnewscale::new_scale_color() 

pt_upgma3 <- pt_upgma2 +
  ggtree::geom_tiplab(
    mapping = aes(color = SpeciesName, label = nodeLabs),
    size = 3
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  scale_color_manual(
    values = c("P. fontis" = "red", "P. brasiliense" = "#E57A44"), breaks = NULL,
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()

ggsave(
  plot = pt_upgma3, filename = file.path(outDir, "upgma_tree.png"),
  width = 10, height = 10
)

#######################################
## function to add annotation to ggtree
annotate_ggtree <- function(pt, offset){
  pt2 <- pt +
    ggtreeExtra::geom_fruit(
      mapping = aes(starshape = type_material),
      geom = "geom_star", fill = "#62BFED", size = 2, starstroke=0.1,
      offset = offset, pwidth = 0.01
    ) +
    ggstar::scale_starshape_manual(
      values = c("type strain" = 1)
    ) + 
    ggnewscale::new_scale_color() +
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "source", color = source),
      geom = "geom_point", shape = 15,
      pwidth = 0.01
    ) +
    scale_color_manual(
      values = c("NCBI" = "#406495", "WUR" = "#468e30", "NVWA" = "#ff6600")
    ) +
    ggnewscale::new_scale_fill() +
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, fill = taxonomy_check_status),
      geom = "geom_star", starshape = 12, size = 2, color = alpha("white", 0),
      pwidth = 0.01
    ) +
    scale_fill_manual(
      values = c("Failed" = "red", "Inconclusive" = "blue", "OK" = alpha("white", 0))
    )
}
#######################################

pt_upgma4 <- annotate_ggtree(pt = pt_upgma3, offset = 0.2) 

# source done
# type_material done
# taxonomy_check_status done
# virulence
# virulence_pcr


ggsave(
  plot = pt_upgma4, filename = file.path(outDir, "upgma_tree.pdf"),
  width = 10, height = 20, scale = 2
)

################################################################################
## highlight outgroup ANI values on ANI score distribution
pt_aniDistr <- dplyr::bind_rows(
  dplyr::slice_sample(.data = aniDf, prop = 0.05),
  dplyr::filter(aniDf, g1 == "104326-106-074" | g1 == "GCF_000803215.1_ASM80321v1")
) %>% 
  dplyr::distinct() %>% 
  ggplot2::ggplot(
    mapping = aes(x = "g1", y = ani, color = g1)
  ) +
  ggbeeswarm::geom_quasirandom(alpha = 0.8) +
  scale_color_manual(
    values = c("104326-106-074" = "blue", "GCF_000803215.1_ASM80321v1" = "red"),
    labels = c(
      "104326-106-074" = sampleInfoList[["104326-106-074"]]$nodeLabs,
      "GCF_000803215.1_ASM80321v1" = sampleInfoList[["GCF_000803215.1_ASM80321v1"]]$nodeLabs
    )
  ) +
  labs(
    title = "highlight outgroup ANI values on ANI score distribution",
    y = "ANI"
  ) +
  pt_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  plot = pt_aniDistr, filename = file.path(outDir, "outgroup_ANI_dist.png"),
  width = 7, height = 7
)

pt_carotovorumLen <- dplyr::filter(sampleInfo, SpeciesName == "P. carotovorum") %>% 
  ggplot(
    mapping = aes(x = "P. carotovorum", y = length)
  ) +
  geom_boxplot(outlier.alpha = 0) +
  ggbeeswarm::geom_quasirandom(mapping = aes(color = id), size = 3) +
  scale_color_manual(
    values = c("104326-106-074" = "blue"),
    labels = c(
      "104326-106-074" = sampleInfoList[["104326-106-074"]]$nodeLabs
    )
  ) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Genome length for all P. carotovorum genomes",
    y = "length"
  ) +
  pt_theme +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    legend.position = "bottom",
    legend.title = element_blank()
  )

ggsave(
  plot = pt_carotovorumLen, filename = file.path(outDir, "P_carotovorum_len.png"),
  width = 7, height = 7
)
################################################################################
## rooting NJ tree with outgroup

outGroup <- dplyr::filter(
  sampleInfo, SpeciesName == "P. fontis"
) %>% 
  dplyr::pull(id)

## negative length edges => 0
treeNj$edge.length[treeNj$edge.length < 0] <- 0

# treedtNj <- as_tibble(treeNj) %>%
#   dplyr::full_join(y = sampleInfo, by = c("label" = "id")) %>%
#   treeio::as.treedata()

treeNjRooted <- root(phy = treeNj, outgroup = outGroup, edgelabel = TRUE)

pt_nj1 <- ggtree::ggtree(tr = treeNjRooted) +
  labs(title = "NJ tree")


pt_nj2 <- pt_nj1 %<+% sampleInfo + 
  geom_point(
    mapping = aes(shape = SpeciesName, color = SpeciesName),
    size = 4
  ) +
  scale_shape_manual(name = "outgroup", values = c("P. fontis" = 16)) +
  scale_color_manual(name = "outgroup", values = c("P. fontis" = "red")) +
  ggnewscale::new_scale_color() 


pt_nj3 <- pt_nj2 +
  ggtree::geom_tiplab(
    mapping = aes(color = SpeciesName, label = nodeLabs),
    size = 3, align = TRUE
  ) +
  scale_color_manual(
    values = c("P. fontis" = "red", "P. brasiliense" = "#E57A44"), breaks = NULL,
    na.value = "black"
  ) +
  ggnewscale::new_scale_color()


pt_nj4 <- annotate_ggtree(pt = pt_nj3, offset = 0.25)


ggsave(
  plot = pt_nj4, filename = file.path(outDir, "nj_tree.pdf"),
  width = 10, height = 20, scale = 2
)


# cophenetic.phylo() computes the pairwise distances between the pairs of tips 
# from a phylogenetic tree using its branch lengths.
# drop.tip() Remove Tips in a Phylogenetic Tree
# node.depth() Depth and Heights of Nodes and Tips
# root()
# ladderize()


## draw tree 

################################################################################










