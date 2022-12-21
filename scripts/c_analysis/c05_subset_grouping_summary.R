suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))


rm(list = ls())

################################################################################

file_fastani <- here::here("analysis", "02_fastANI", "ANI_results")
file_metadata <- here::here("data/pangenomes/pectobacterium.v2", "genomes_metadata.csv")
file_sampleSummary <- here::here("data/reference_data", "sample_metadata.tsv")
file_tree <- here::here("analysis", "02_fastANI", "ANI_UPGMA.newick")

outDir <- here::here("analysis", "03_pangenome_pecto_v2", "subset_optimal_group")

files_grpOverview <- list.files(
  path = "analysis/03_pangenome_pecto_v2/subset_optimal_group",
  pattern = "fn_fp.csv", recursive = TRUE, full.names = TRUE
)

################################################################################
df <- suppressMessages(
  readr::read_csv(
    file = "analysis/03_pangenome_pecto_v2/subset_optimal_group/optimal_grouping_pangenome/grouping_overview.csv"
  )
)

fpFn <- purrr::map_dfr(
  .x = files_grpOverview,
  .f = function(fl){
    suppressMessages(readr::read_csv(file = fl)) %>% 
      dplyr::mutate(
        id = stringr::str_replace(
          string = fl, 
          pattern = ".*optimal_grouping_(\\w+)/fn_fp.csv",
          replacement = "\\1"
        )
      )
  }
) %>% 
  dplyr::mutate(
    setType = dplyr::case_when(
      id == "pangenome" ~ "pangenome",
      id == "typeStrain" ~ "type_strains",
      TRUE ~ "random"
    ),
    setType = forcats::fct_relevel(.f = setType, "pangenome", "type_strains")
  )


pt_optGrp <- ggplot(
  # data = dplyr::filter(fpFn, id == "rand_004"),
  data = fpFn,
  mapping = aes(x=Mode, y=Score, group = id)
) +
  geom_line(
    mapping = aes(color=Type, linetype = setType, alpha = setType, size = setType)
  ) +
  scale_color_manual(values=c("FP" = "#000000", "FN" = "#990000")) +
  scale_linetype_manual(
    values = c("random" = "dotted", "type_strains" = "solid", "pangenome" = "solid")
  ) +
  scale_alpha_manual(values = c("random" = 0.7, "type_strains" = 0.7, "pangenome" = 1)) +
  scale_size_manual(values = c("random" = 0.5, "type_strains" = 1, "pangenome" = 2)) +
  scale_y_continuous(
    name = "Proteins", trans = 'log10',
    breaks=c(1,10,100,1000,10000,100000,100000,1000000)
  ) +
  scale_x_continuous("Relaxation mode", breaks = 1:8, limits=c(1, 8)) +
  labs(
    title = "Optimized grouping with subsets"
  ) +
  facet_grid(rows = vars(Type)) +
  theme_bw(base_size = 15) +
  theme(
    legend.key.size = unit(1.5, "cm"),
    legend.position = "right"
  )


ggsave(
  plot = pt_optGrp, width = 10, height = 8,
  filename = file.path(outDir, "subset_optimal_groups.png")
)








