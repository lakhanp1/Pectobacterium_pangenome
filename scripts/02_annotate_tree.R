suppressMessages(library(ggtree))
suppressMessages(library(treeio))
suppressMessages(library(tidyverse))
suppressMessages(library(RColorBrewer))
suppressMessages(library(viridis))
suppressMessages(library(ggpubr))
suppressMessages(library(here))
suppressMessages(library(ComplexHeatmap))

## plot tree and phenotype annotations

rm(list = ls())

##################################################################################
analysisName <- "ANI_tree"
outDir <- here::here("analysis", "04_pangenome.50g")

outPrefix <- file.path(outDir, analysisName)

pangenome <- here::here("analysis", "04_pangenome.50g", "pectobacterium.50g.DB")

file_phenotypes <- file.path(outDir, "pectobacterium.50g.phenotypes.txt")
file_aniTree <- file.path(pangenome, "ANI", "fastANI", "ANI_RENAMED.newick")
file_aniScores <- file.path(pangenome, "ANI", "fastANI", "ANI_scores.csv")

##################################################################################

phenotypes <- suppressMessages(readr::read_tsv(file_phenotypes)) %>% 
  tidyr::unite(genomeId, Genome, id, sep = "_", remove = F) %>% 
  tidyr::unite(pgSampleId, Genome, species, sep = " ", remove = F) %>% 
  tidyr::unite(pgSampleId, pgSampleId, batch, sep = "_", remove = F) %>% 
  dplyr::mutate(
    pgSampleId = stringr::str_replace(
      string = pgSampleId, pattern = "_old", replacement = ""
    )
  )

aniTreeTibble <- treeio::read.newick(file = file_aniTree) %>% 
  # treeio::as.treedata() %>% 
  as_tibble()

aniScores <- suppressMessages(read_csv(file = file_aniScores)) %>% 
  dplyr::left_join(
    y = dplyr::select(phenotypes, genomeId, id, pgSampleId),
    by = c("Genome" = "pgSampleId")
  ) %>% 
  dplyr::select(genomeId, id, everything(), -Genome) %>% 
  dplyr::rename(!!! set_names(x = phenotypes$pgSampleId, nm = phenotypes$genomeId))
  

aniTree <- dplyr::full_join(
  x = aniTreeTibble, y = phenotypes, by = c("label" = "genomeId")
) %>% 
  treeio::as.treedata()

# as_tibble(aniTree) %>% view()
# 
# treeio::get.fields(aniTree)
# treeio::get.data(aniTree) %>% view()
##################################################################################



pt_tree <- ggtree::ggtree(
  tr = aniTree,
  # mapping = aes(color = species)
  # branch.length = "none"
) +
  scale_x_continuous(
    expand = expansion(add = c(0.01, 0.06))
  ) +
  geom_tiplab(
    mapping = aes(color = species, label = id),
    align = TRUE
  ) +
  geom_treescale(x = 0, y = length(aniTree@phylo$tip.label)-2) 

leafOrder <- dplyr::arrange(.data = pt_tree$data, y) %>% 
  dplyr::filter(isTip) %>% 
  dplyr::pull(label)


## IMP: set the factor levels so rest of the plots have same Y axis order as tree leaves
phenotypes <- dplyr::mutate(
  .data = phenotypes,
  genomeId = forcats::fct_relevel(genomeId, !!!leafOrder),
  pcr = dplyr::case_when(
    pcr == "positive" ~ "+",
    pcr == "negative" ~ "-",
    TRUE ~ pcr
  )
)

##################################################################################
## ANI score heatmap with tree

tempTree <- treeio::read.newick(file = file_aniTree)
ape::is.ultrametric(tempTree)
# ape::chronos(tempTree)
# as.hclust(tempTree)

aniScoresMat <- dplyr::select(aniScores, -id)%>% 
  tibble::column_to_rownames(var = "genomeId") %>% 
  as.matrix()

pt_aniHt <- ComplexHeatmap::Heatmap(
  matrix = aniScoresMat,
  col = circlize::colorRamp2(
    breaks = seq(from = 91, length.out = 9, by = 1),
    colors = RColorBrewer::brewer.pal(n = 9, name = "RdPu")
  ) ,
  cluster_rows = F,
  row_order = rev(leafOrder),
  column_order = rev(leafOrder),
  row_labels = aniScores$id,
  column_title = "Pairwise ANI comparison",
  column_title_gp = gpar(fontsize = 20, fontface = "bold"),
  heatmap_legend_param = list(title = "ANI %")
)

png(filename = paste(outPrefix, ".score_heatmap.png", sep = ""), width = 1500, height = 1500, res = 150)
pt_aniHt
dev.off()

##################################################################################
## generate annotation plots
pt_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(size = 14, angle = 90, face = "bold", vjust = 0.5),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank()
  )

pt_symptoms <- dplyr::select(phenotypes, genomeId, starts_with("symptom_per")) %>% 
  tidyr::pivot_longer(
    cols = -genomeId
  ) %>% 
  dplyr::mutate(
    name = stringr::str_replace(string = name, pattern = "symptom_per_", replacement = "")
  ) %>% 
  ggplot2::ggplot(
    mapping = aes(x = name, y = genomeId, fill = value)
  ) +
  geom_tile() +
  labs(title = "Symptom\n%") +
  scale_fill_viridis(
    option="D", na.value = "grey80",
    name = "Symptom %"
  ) +
  pt_theme +
  theme(
    plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0)
  )

pt_vir <- dplyr::select(phenotypes, genomeId, virulence) %>% 
  tidyr::pivot_longer(
    cols = -genomeId
  ) %>% 
  dplyr::mutate(
    value = forcats::fct_relevel(
      value, "yes", "maybe", "no", "unknown"
    )
  ) %>% 
  ggplot2::ggplot(
    mapping = aes(x = name, y = genomeId, shape = value, color = value)
  ) +
  geom_point(size = 3, stroke = 2) +
  # labs(title = "Virulence") +
  scale_shape_manual(
    name = "Virulence",
    values = c(
      "yes" = 17, "maybe" = 2, "no" = 19, "unknown" = 4
    )
  ) +
  scale_color_manual(
    name = "Virulence",
    values = c(
      "yes" = "red", "maybe" = "red", "no" = "green", "unknown" = "black"
    )
  ) +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )


pt_pcr <- dplyr::select(phenotypes, genomeId, pcr) %>% 
  tidyr::pivot_longer(
    cols = -genomeId
  ) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "PCR", y = genomeId, label = value, color = value)
  ) +
  geom_text(size = 8, fontface = "bold") +
  scale_color_manual(
    values = c("+" = "red", "-" = "green"), breaks = NULL
  ) +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )

pt_sampleYear <- dplyr::select(phenotypes, genomeId, sampling_year) %>% 
  ggplot2::ggplot(
    mapping = aes(x = sampling_year, y = genomeId)
  ) +
  geom_point(size = 4) +
  scale_x_continuous(expand = expansion(add = 0.5)) +
  labs(title = "Sampling\nyear") +
  pt_theme +
  theme(
    plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0)
  )

pt_tissue <- dplyr::select(phenotypes, genomeId, tissue) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "tissue", y = genomeId, color = tissue)
  ) +
  geom_point(size = 4, shape = 15) +
  scale_color_manual(
    values = c("Stem" = "seagreen4", "Tuber" = "wheat3"),
    breaks = c("Stem", "Tuber"), name = "Tissue"
  ) +
  labs(title = "Tissue") +
  pt_theme +
  theme(
    plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )

pt_region <- dplyr::select(phenotypes, genomeId, region) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "Region", y = genomeId, label = region)
  ) +
  geom_text(size = 4) +
  labs(title = "Region") +
  pt_theme +
  theme(
    plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()
  )

## arrange plots one by one
pt_all <- pt_vir %>% aplot::insert_left(pt_tree, width = 14) %>% 
  aplot::insert_right(pt_symptoms, width = 2) %>% 
  aplot::insert_right(pt_pcr, width = 1) %>% 
  aplot::insert_right(pt_sampleYear, width = 2) %>% 
  aplot::insert_right(pt_tissue, width = 1.5) %>% 
  aplot::insert_right(pt_region, width = 3)


## arrange plots as a list
pt_all2 <- aplot::plot_list(
  gglist = list(pt_tree, pt_vir, pt_symptoms, pt_pcr, pt_sampleYear, pt_tissue, pt_region),
  nrow = 1,
  widths = c(15, 1, 4, 1, 6, 1, 3),
  guides = "collect"
)

png(filename = paste(outPrefix, ".png", sep = ""), width = 2500, height = 1500, res = 160)
pt_all
dev.off()

