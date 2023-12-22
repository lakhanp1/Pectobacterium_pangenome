#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))

# extract homology groups for plasmid regions

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

################################################################################
sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId 
)

plasmidDf <- suppressMessages(
  readr::read_tsv(confs$data$plasmids$files$data)
)

################################################################################
# panPlasmids <- dplyr::left_join(sampleInfo, plasmidDf, by = "sampleId")

dplyr::group_by(plasmidDf, sampleId, SpeciesName) %>%
  dplyr::summarise(n = n_distinct(plasmid_id, na.rm = TRUE), .groups = "drop") %>%
  dplyr::group_by(SpeciesName) %>%
  skimr::skim()

# get homology groups for each plasmid region
plasmidHgs <- dplyr::select(plasmidDf, sampleId, genomeId, chr, plasmid_id) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    hgs = list(
      region_homology_groups(
        pandb = orgDb, genome = genomeId, chr = chr
      )
    )
  )

# filter to remove prophages without any hgs 
plasmidHgs %>%
  dplyr::mutate(
    nHgs = length(hgs),
    hgs = paste(hgs, collapse = ";")
  ) %>%
  dplyr::filter(nHgs > 0) %>%
  dplyr::select(plasmid_id, sampleId, nHgs, hgs) %>%
  readr::write_tsv(
    file = confs$analysis$plasmids$preprocessing$files$raw_plasmid_hg
  )

################################################################################

plasmidHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$plasmids$preprocessing$files$raw_plasmid_hg)
) %>% 
  dplyr::left_join(
    y = dplyr::select(sampleInfo, sampleId, genomeId, SpeciesName),
    by = "sampleId"
  )

## binary matrix for homology_group PAV
hgBinaryMat <- homology_groups_mat(pandb = panOrgDb, type = "pav")

sppGrpStats <- NULL

## get species wise core, accessory, unique group stats and GO
for (sp in c(unique(sampleInfo$SpeciesName), "Pangenome")) {
  
  atPangenomeScale <- sp == "Pangenome"
  
  if (atPangenomeScale) {
    spGenomes <- sampleInfo$genomeId
  } else{
    spGenomes <- dplyr::filter(sampleInfo, SpeciesName == .env$sp) %>%
      dplyr::pull(genomeId)
  }
  
  # prophage homology groups: handle special case when HGs for Pangenome are needed
  spRegionHgs <- dplyr::filter(
    plasmidHgs, SpeciesName == !!sp | atPangenomeScale
  ) %>% 
    dplyr::select(hgId = hgs) %>% 
    dplyr::mutate(hgId = stringr::str_split(hgId, ";")) %>% 
    tidyr::unnest(cols = c(hgId)) %>% 
    dplyr::distinct() %>%
    dplyr::mutate(
      hgInRegion = "plasmid",
      SpeciesName = .env$sp
    )
  
  # pangenome homology groups
  spHgTypes <- matrixStats::colSums2(
    x = hgBinaryMat, useNames = T,
    rows = which(rownames(hgBinaryMat) %in% spGenomes)
  ) %>%
    tibble::enframe(name = "hgId", value = "nGenomes") %>%
    dplyr::filter(nGenomes != 0) %>%
    dplyr::mutate(
      class = dplyr::case_when(
        nGenomes == 1 ~ "unique",
        nGenomes == !!length(spGenomes) ~ "core",
        nGenomes < !!length(spGenomes) & nGenomes > 1 ~ "accessory"
      )
    ) %>% 
    dplyr::left_join(y = spRegionHgs, by = "hgId") %>% 
    tidyr::replace_na(list(hgInRegion = "non-plasmid", SpeciesName = sp))
  
  # if only one genome for species, all HGs should be assigned to "core"
  if(length(spGenomes) == 1){
    spHgTypes$class <- "core"
  }
  
  # phageHgDf <- dplyr::filter(spHgTypes, hgInRegion == "plasmid") %>% 
  #   dplyr::select(hgId, !!sp := class) %>% 
  #   dplyr::full_join(phageHgDf, by = "hgId")
  
  ## group stats
  proHgSpStats <- dplyr::count(spHgTypes, class, hgInRegion, name = "plasmid_nHgs") %>%
    dplyr::filter(hgInRegion == "plasmid") %>% 
    dplyr::left_join(
      y = as.data.frame(table(spHgTypes$class), responseName = "nHgs"),
      by = c("class" = "Var1")
    ) %>% 
    dplyr::select(-hgInRegion) %>% 
    dplyr::bind_rows(
      dplyr::summarise(
        .,
        dplyr::across(.cols = -class, .fns = sum),
        dplyr::across(.cols = class, .fns = ~'total')
      )
    ) %>% 
    dplyr::mutate(
      SpeciesName = .env$sp,
      plasmidRatio = round(plasmid_nHgs / nHgs, digits = 3)
    ) %>%
    dplyr::select(SpeciesName, class, nHgs, plasmid_nHgs, plasmidRatio)
  
  sppGrpStats <- dplyr::bind_rows(sppGrpStats, proHgSpStats)
}

spec1 <- tidyr::build_wider_spec(
  data = sppGrpStats,
  names_from = class,
  values_from = c(nHgs, plasmid_nHgs, plasmidRatio),
  names_glue = "{.value}.{class}"
) %>% 
  dplyr::mutate(
    .name = dplyr::if_else(
      .value == "nHgs", true = class, false = .name
    )
  )

plasmidSpeciesStats <- tidyr::pivot_wider_spec(data = sppGrpStats, spec = spec1)

plasmidSpeciesStats %<>% dplyr::select(
  SpeciesName, core, accessory, unique, total,
  ends_with(c(".total", ".core", ".accessory", ".unique"))
) %>% 
  dplyr::arrange(desc(total))

plasmidSpeciesStats %>% 
  readr::write_tsv(
    file = confs$analysis$plasmids$files$stats_species
  )
