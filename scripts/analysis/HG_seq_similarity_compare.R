#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(optparse))

## get pairwise sequence similarity comparison of two genomes for two sets of
## homology groups:
## 1) group specific for a phenotype
## 2) other homology groups

## This script will be replaced once I develop API to access information from pangenome
## and store it in efficient data structure for local analysis

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser,
  opt_str = c("-1", "--g1"), type = "character", action = "store",
  help = "genome 1 ID"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-2", "--g2"), type = "character", action = "store",
  help = "comma separated list of genome IDs to be compared against genome 1"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("--phenotype"), type = "character", action = "store",
  help = "phenotype"
)

parser <- optparse::add_option(
  parser,
  opt_str = c("-c", "--config"), type = "character", action = "store",
  help = "project config YAML file"
)

opts <- optparse::parse_args(parser)

if (any(is.null(opts$g1), is.null(opts$g2), is.null(opts$config))) {
  stop(optparse::print_help(parser), call. = TRUE)
}

# ## for test
# opts$config <- "project_config.yaml"
# opts$g1 <- "399"
# opts$g2 <- "404,405,406,18,26,323,334,93,422,10,9,211,304,285,437,445,33,201,224,276,312"
# opts$phenotype <- "assay_FN"
# #######

################################################################################
confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = opts$config)),
  dir = "."
)

g1 <- opts$g1
g2 <- unlist(stringr::str_split(string = opts$g2, pattern = ","))

phenotype <- opts$phenotype

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]


outDir <- file.path(confs$analysis$association$path, phenotype)
outPrefix <- paste(outDir, "/", phenotype, ".hg_seq_identity.", g1, sep = "")

## sequence info file for mRNAs across all genomes belonging to homology groups
## specific to a phenotype of interest
# file_phenoSpecificHg <- confs$analysis$association$files$pheno_specific_groups

################################################################################
#####***********
#####* in future, remove this file dependency
#####**********
phenoSpecificGrp <- suppressMessages(
  readr::read_tsv(
    file = confs$analysis$association$files$pheno_specific_groups,
    col_types = "cc", col_names = c("phenotype", "group_id")
  )
) %>%
  dplyr::filter(phenotype == !!phenotype) %>%
  dplyr::mutate(group_id = stringr::str_split(group_id, ",")) %>%
  tidyr::unnest(cols = c(group_id))


## import homology group data
hgData <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups
  )
) %>%
  dplyr::rename_with(.fn = ~ stringr::str_replace_all(.x, "( |-)", "_")) %>%
  dplyr::rename_with(.fn = ~ stringr::str_replace_all(.x, "Genome_", "")) %>%
  dplyr::rename_with(.fn = ~ tolower(.x)) %>%
  dplyr::select(group_id = homology_group_id, class, all_of(c(g1, g2))) %>%
  dplyr::mutate(group_id = as.character(group_id))

## get homology groups shared by the two genomes
filteredHgs <- dplyr::filter(
  hgData,
  !!sym(g1) != 0,
  dplyr::if_any(.cols = !!g2, .fns = ~ . != 0)
) %>%
  dplyr::left_join(phenoSpecificGrp, by = "group_id") %>%
  dplyr::mutate(
    type = dplyr::if_else(is.na(phenotype), "other", phenotype),
    class = stringr::str_replace(class, "single copy orthologous", "SCO")
  )

# sharedHgs <- dplyr::bind_rows(
#   dplyr::filter(filteredHgs, type == "other") %>%
#     dplyr::slice_sample(n = 500),
#   dplyr::filter(filteredHgs, type != "other")
# )


## pairwise genome similarity
groupSeqSimilarity <- tibble::tibble(
  group_id = character(), g1 = character(), g2 = character(),
  g1mRNA = character(), g2mRNA = character(), identity = double(),
  class = character(), type = character()
)

readr::write_tsv(
  groupSeqSimilarity,
  file = paste(outPrefix, ".data.tab", sep = "")
)

# hg <- "22573071"
## get mRNA-id
for (hg in filteredHgs$group_id) {
  cat(hg, " ")

  file_grpSeqInfo <- file.path(
    panConf$db$path, "alignments/msa_per_group/grouping_v4", hg, "input", "sequences.info"
  )

  if (!file.exists(file_grpSeqInfo)) {
    next
  }

  mRNAs <- suppressMessages(readr::read_csv(file_grpSeqInfo, skip = 5)) %>%
    dplyr::select(
      Genome = "#genome", mRNA_id = "mRNA identifier"
    ) %>%
    dplyr::filter(Genome %in% c(g1, g2)) %>%
    dplyr::pull(mRNA_id, name = "Genome")

  # safety
  if (length(mRNAs) < 2) {
    next
  }

  file_similarity <- file.path(
    panConf$db$path, "/alignments/msa_per_group/grouping_v4", hg,
    "output/similarity_identity", "nuc_trimmed_identity.csv"
  )

  # safety
  if (!file.exists(file_similarity)) {
    next
  }

  g1mRNA <- unname(mRNAs[g1])

  perIdentity <- suppressMessages(readr::read_csv(file_similarity)) %>%
    dplyr::select(Sequences, !!unname(mRNAs[g1])) %>%
    dplyr::rename(g2mRNA = Sequences, identity = !!unname(mRNAs[g1]))

  hgSimilarity <- tibble::tibble(
    group_id = hg,
    g1 = g1, g2 = g2, g1mRNA = unname(mRNAs[g1]), g2mRNA = unname(mRNAs[g2])
  ) %>%
    dplyr::left_join(y = perIdentity, by = "g2mRNA") %>%
    dplyr::bind_rows(groupSeqSimilarity) %>%
    dplyr::left_join(
      y = dplyr::select(filteredHgs, group_id, class, type), by = "group_id"
    )

  readr::write_tsv(
    hgSimilarity,
    file = paste(outPrefix, ".data.tab", sep = ""),
    col_names = FALSE, append = TRUE
  )
}

################################################################################
# plotting the data
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

tileOrder <- dplyr::filter(filteredHgs, type == !!phenotype) %>%
  dplyr::select(all_of(g2)) %>%
  dplyr::summarise(dplyr::across(.fns = ~ sum(.x))) %>%
  unlist() %>%
  sort(decreasing = TRUE) %>%
  names()

hgSimilarity <- suppressMessages(
  readr::read_tsv(file = paste(outPrefix, ".data.tab", sep = ""), col_types = "ccc")
) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, g2Species = SpeciesName),
    by = c("g2" = "Genome")
  ) %>%
  dplyr::mutate(
    class = forcats::fct_relevel(class, "accessory", "core & SCO", "core"),
    g2 = forcats::fct_relevel(g2, !!!tileOrder)
  )

pt_identity <- ggplot(
  hgSimilarity,
  mapping = aes(x = type, y = identity)
) +
  ggbeeswarm::geom_quasirandom(alpha = 0.7, color = "black") +
  coord_cartesian(ylim = c(50, 100)) +
  labs(
    title = "homology group pairwise similarity",
    subtitle = paste(
      "pairwise similarity of homology groups common between",
      "genome", g1, "(", sampleInfo$SpeciesName[which(sampleInfo$Genome == g1)],
      ") and other genomes"
    ),
    y = "%identity"
  ) +
  # scale_color_viridis_d() +
  facet_wrap(
    facets = ~ g2 + g2Species,
    labeller = label_wrap_gen(multi_line = FALSE)
  ) +
  theme_bw(base_size = 12) +
  theme(
    # legend.position = "bottom",
    panel.grid.major.x = element_blank()
  )

ggsave(
  filename = paste(outPrefix, ".png", sep = ""), plot = pt_identity, width = 10, height = 8
)
