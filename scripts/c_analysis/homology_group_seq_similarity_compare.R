#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(optparse))

## plot pairwise sequence similarity comparison of two genomes for two sets of 
## homology groups:
## 1) group specific for a phenotype
## 2) other homology groups

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/01_RScripts/02_R_utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

## argment parsing
parser <- optparse::OptionParser()

parser <- optparse::add_option(
  parser, opt_str = c("-1", "--g1"), type = "character", action = "store",
  help = "genome 1 ID"
)

parser <- optparse::add_option(
  parser, opt_str = c("-2", "--g2"), type = "character", action = "store",
  help = "genome 1 ID"
)

parser <- optparse::add_option(
  parser, opt_str = c("-c", "--config"), type = "character", action = "store",
  help = "project config YAML file"
)

opts <- optparse::parse_args(parser)

if(any(is.null(opts$g1), is.null(opts$g2), is.null(opts$config)))
  stop(optparse::print_help(parser), call. = TRUE)

# ## for test
# opts$config <- "project_config.yaml"
# opts$g1 <- "399"
# opts$g2 <- "18"
# #######

################################################################################
confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = opts$config)),
  dir = "."
)

g1 <- opts$g1
g2 <- opts$g2
phenotype <- "assay_FN"
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]


outDir <- file.path(confs$analysis$association$dir, phenotype)
outPrefix <- paste(outDir, "/", phenotype, ".hg_seq_identity.", g1, "_", g2, sep = "")

## sequence info file for mRNAs across all genomes belonging to homology groups 
## specific to a phenotype of interest
file_seqInfo <- paste(outPrefix, ".pheno_specific.seq_info.txt", sep = "")
# file_phenoSpecificHg <- confs$analysis$association$files$pheno_specific_groups

################################################################################

phenoSpecificGrp <- suppressMessages(
  readr::read_tsv(
    confs$analysis$association$files$pheno_specific_groups,
    col_types = "cc", col_names = c("phenotype", "group_id"))
) %>% 
  dplyr::filter(phenotype == !!phenotype) %>% 
  dplyr::mutate(group_id = stringr::str_split(group_id, ",")) %>% 
  tidyr::unnest(cols = c(group_id))


## import homology group data
g2hg <- suppressMessages(
  readr::read_csv(
    file = panConf$db$gene_classification$GC.100.0$files$groups
  )
) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
  dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "Genome_", "")) %>% 
  dplyr::rename_with(.fn = ~tolower(.x)) %>% 
  dplyr::select(group_id = homology_group_id, class, all_of(c(g1, g2))) %>% 
  dplyr::mutate(group_id = as.character(group_id))

## get homology groups shared by the two genomes
filteredHgs <- dplyr::filter(
  g2hg,
  !!sym(g1) != 0, !!sym(g2) != 0
) %>% 
  dplyr::left_join(phenoSpecificGrp, by = "group_id") %>% 
  dplyr::mutate(
    type = dplyr::if_else(is.na(phenotype), "other", phenotype),
    class = stringr::str_replace(class, "single copy orthologous", "SCO")
  )

sharedHgs <- dplyr::bind_rows(
  dplyr::filter(filteredHgs, type == "other") %>% 
    dplyr::slice_sample(n = 250),
  dplyr::filter(filteredHgs, type != "other")
)


## pairwise genome similarity
groupSeqSimilarity <- tibble::tibble()

# hg <- "22579048"
## get mRNA-id
for (hg in sharedHgs$group_id) {
  
  cat(hg, " ")
  
  file_grpSeqInfo <- file.path(
    panConf$db$dir, "alignments/msa_per_group/grouping_v4", hg, "input", "sequences.info"
  )
  
  if(!file.exists(file_grpSeqInfo)){
    next
  }
  
  mRNAs <- suppressMessages(readr::read_csv(file_grpSeqInfo, skip = 5)) %>% 
    dplyr::select(
      Genome = "#genome", mRNA_id = "mRNA identifier"
    ) %>% 
    dplyr::filter(Genome %in% c(g1, g2)) %>% 
    dplyr::pull(mRNA_id, name = "Genome")
  
  # safety
  if(length(mRNAs) != 2){
    next
  }
  
  file_similarity <- file.path(
    panConf$db$dir, "/alignments/msa_per_group/grouping_v4", hg,
    "output/similarity_identity", "nuc_trimmed_identity.csv"
  )
  
  # safety
  if(!file.exists(file_similarity)){
    next
  }
  
  perIdentity <- suppressMessages(readr::read_csv(file_similarity)) %>% 
    dplyr::select(Sequences, !!unname(mRNAs)) %>% 
    dplyr::filter(Sequences %in% !!mRNAs) %>% 
    tibble::column_to_rownames(var = "Sequences") %>% 
    as.matrix()
  
  groupSeqSimilarity <- tibble::tibble(
    group_id = hg,
    g1 = g1, g2 = g2, g1mRNA = unname(mRNAs[g1]), g2mRNA = unname(mRNAs[g2]),
    identity = perIdentity[unname(mRNAs[g1]), g2mRNA = unname(mRNAs[g2])]
  ) %>% 
    dplyr::bind_rows(groupSeqSimilarity)
  
}


sharedHgs %<>% dplyr::left_join(groupSeqSimilarity, by = "group_id") %>% 
  dplyr::filter(!is.na(identity))

readr::write_tsv(sharedHgs, file = paste(outPrefix, ".data.tab", sep = ""))

pt_identity <- ggplot(
  sharedHgs,
  mapping = aes(x = type, y = identity, color = class)
) +
  ggbeeswarm::geom_quasirandom() +
  # geom_jitter() +
  labs(
    title = "homology group pairwise similarity",
    y = "%identity"
  ) +
  scale_color_viridis_d() +
  theme_bw(base_size = 16) 

ggsave(
  filename = paste(outPrefix, ".png", sep = ""), plot = pt_identity, width = 6, height = 6
)











