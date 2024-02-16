suppressMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressMessages(library(tidyverse))

# Extract homology groups for one genome
# Use nodes on phylogenetic tree to select a set of genomes based on which
# the core, accessory and unique categories will be assigned.

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")

################################################################################
genomeOfInterest <- "g_182"

treeMethod <- "kmer_upgma" # ani_upgma, kmer_upgma
nodeInclude <- "n11"
nodeExclude <- "n145"

panOrgDb <- org.Pectobacterium.spp.pan.eg.db

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

column_nodepath <- paste("nodepath.", treeMethod, sep = "")
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

subGenomes <- dplyr::filter(
  sampleInfo,
  stringr::str_detect(
    string = !!sym(column_nodepath), pattern = paste(";", nodeInclude, ";", sep = "")
  ),
  stringr::str_detect(
    string = !!sym(column_nodepath), pattern = paste(";", nodeExclude, ";", sep = ""),
    negate = TRUE
  )
)

################################################################################
# extract data
gn <- AnnotationDbi::select(
  x = panOrgDb, keys = genomeOfInterest, keytype = "genomeId",
  columns = c(
    "GID", "class", "chr_name", "start", "end", "strand", "chr", "mRNA_id",
    "genePos", "COG_description", "pfam_description"
  )
) %>%
  dplyr::mutate(
    dplyr::across(.cols = c(start, end), .fns = as.numeric),
  ) %>%
  dplyr::group_by(mRNA_id, GID) %>%
  dplyr::mutate(
    dplyr::across(
      .cols = c(COG_description, pfam_description),
      .fns = ~ paste(unique(.x), collapse = "; ")
    )
  ) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(chr_name, start) %>% 
  tidyr::unite(
    col = cog_pfam, COG_description, pfam_description, sep = "; ", remove = FALSE
  ) %>% 
  dplyr::mutate(
    broad_function = stringr::str_extract_all(
      string = cog_pfam,
      pattern = regex(
        "integrase|secretion.system|methyltransferase|restriction|integrase|transposase|toxin",
        ignore_case = TRUE
        )
    ),
    broad_function = purrr::map_chr(
      .x = broad_function,
      .f = function(x){
        if(length(x) > 0){
          return(paste(unique(tolower(x)), collapse = ";"))
        } else {
          return("")
        }
      }
    )
  ) %>% 
  dplyr::select(-cog_pfam)

# HG PAVs for sub pangenome
subPanHgs <- sub_pangenome_hgs(
  pandb = panOrgDb,
  genomes = union(genomeOfInterest, subGenomes$genomeId)
) %>% 
  dplyr::rename(subset_class = class, subset_nGenomes = nGenomes)

# HG PAVs for the pangenome
panHgs <- sub_pangenome_hgs(
  pandb = panOrgDb, genomes = sampleInfo$genomeId
) %>% 
  dplyr::rename(pan_class = class, pan_nGenomes = nGenomes)

# combine data
df <- dplyr::left_join(x = gn, y = subPanHgs, by = c("GID" = "hgId")) %>% 
  dplyr::left_join(y = panHgs, by = c("GID" = "hgId")) %>% 
  dplyr::select(
    hgId = GID, class, pan_class, pan_nGenomes, subset_class, subset_nGenomes, 
    chr_name, start, end, strand, mRNA_id, genePos,
    broad_function, COG_description, pfam_description
  )

clipr::write_clip(df)

################################################################################

