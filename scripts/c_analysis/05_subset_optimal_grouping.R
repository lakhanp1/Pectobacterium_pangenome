suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(ggnewscale))
suppressPackageStartupMessages(library(ape))

## perform optimal grouping on subset of the data:
## 1) 20 clades of UPGMA tree
## 2) only type strains

rm(list = ls())

################################################################################
set.seed(124)


file_fastani <- here::here("analysis", "02_fastANI", "ANI_results")
file_metadata <- here::here("data/pangenomes/pectobacterium.v2", "genomes_metadata.csv")
file_sampleSummary <- here::here("data/reference_data", "sample_metadata.tsv")
file_tree <- here::here("analysis", "02_fastANI", "ANI_UPGMA.newick")

outDir <- here::here("analysis", "03_pangenome_pecto_v2", "subset_optimal_group")

################################################################################

if(!dir.exists(outDir)){
  dir.create(outDir, recursive = T)
}

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

treeUpgma <- ape::read.tree(file = file_tree)
################################################################################

typeStrains <- dplyr::filter(sampleInfo, !is.na(type_material)) %>% 
  dplyr::group_by(SpeciesName) %>% 
  dplyr::arrange(desc(N50)) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()


################################################################################

treeSets <- tibble::enframe(
  x = cutree(tree = as.hclust.phylo(treeUpgma), k = 20),
  name = "id", value = "group"
) %>% 
  dplyr::add_count(group, name = "count") %>% 
  dplyr::arrange(count) %>% 
  dplyr::group_by(group) %>% 
  dplyr::slice_sample(n = 100, replace = TRUE) %>% 
  dplyr::mutate(
    set = 1:n(),
    setId = paste("rand", sprintf(fmt = "%03d", set), sep = "_")
  ) %>% 
  dplyr::ungroup() %>% 
  dplyr::left_join(y = dplyr::select(sampleInfo, id, Genome), by = "id")


optGrpSets <- dplyr::group_by(treeSets, setId) %>% 
  dplyr::summarise(genomes = paste(Genome, collapse = ",")) %>% 
  dplyr::bind_rows(
    tibble(
      setId = "typeStrain", genomes = paste(typeStrains$Genome, collapse = ",")
    )
  ) %>% 
  dplyr::mutate(
    setId = forcats::fct_relevel(setId, "typeStrain")
  ) %>% 
  dplyr::arrange(setId)

readr::write_tsv(
  x = optGrpSets, col_names = F,
  file = file.path(outDir, "subsets_conf.tab")
)


