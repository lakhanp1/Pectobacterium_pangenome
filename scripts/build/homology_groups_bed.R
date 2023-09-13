suppressMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressMessages(library(rtracklayer))

## prepare a homology groups BED file to visualize in IGV

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)
orgDb <- org.Pectobacterium.spp.pan.eg.db
genomeId <- "g_182"

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId 
)

outPrefix <- sampleInfo$sampleId[which(sampleInfo$genomeId == genomeId)]
################################################################################
# extract data
gn <- AnnotationDbi::select(
  x = orgDb, keys = genomeId, keytype = "genomeId",
  columns = c(
    "GID", "class", "chr_name", "start", "end", "strand", "chr", "mRNA_id",
    "COG_description"
  )
) %>%
  dplyr::group_by(mRNA_id, GID) %>%
  dplyr::mutate(
    COG_description = paste(COG_description, collapse = "; ")
  ) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup()

# make GRanges object
gr <- dplyr::mutate(
  gn,
  score = dplyr::case_when(
    class == "core" ~ 1000,
    class == "core & single copy orthologous" ~ 1000,
    class == "accessory" ~ 700,
    class == "unique" ~ 300
  ),
  thickStart = as.integer(start),
  thickEnd = as.integer(end),
  itemRgb = dplyr::case_when(
    class == "core" ~ "#009E73",
    class == "core & single copy orthologous" ~ "#009E73",
    class == "accessory" ~ "#0072B2",
    class == "unique" ~ "#D55E00"
  ),
  name = stringr::str_c(
    "Homology group=", GID, ";COG=", COG_description,
    sep = ""
  ),
  name = stringr::str_replace_all(
    string = name, pattern = "\\s+", replacement = "%20"
  ),
  name = stringr::str_replace_all(
    string = name, pattern = ":", replacement = "-"
  )
) %>%
  dplyr::select(
    chr = chr_name, start, end, strand, score, name,
    thickStart, thickEnd, itemRgb
  ) %>%
  GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE)


gru <- as(gr, "UCSCData")
gru@trackLine <- as(
  "track name=\"homology_groups\" visibility=2 itemRgb=\"On\"",
  "BasicTrackLine"
)

rtracklayer::export.bed(object = gru, con = paste(outPrefix, ".hg.bed", sep = ""))
