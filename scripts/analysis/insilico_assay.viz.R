suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(GenomicRanges))

## parse PCR primer blastn results and generate summary figure

rm(list = ls())


source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
outGroup <- confs$analysis$phylogeny$outgroup

analysisName <- "bioassay_probe_match"

outDir <- confs$analysis$insilico_assay$dir
outPrefix <- file.path(outDir, analysisName)

################################################################################

genomeChrs <- suppressMessages(
  readr::read_tsv(
    file = confs$data$pangenomes[[pangenome]]$files$genome_chrs,
    col_names = c("Genome", "chr"), col_types = "c"
  )
)

sampleInfo <- get_metadata(file = confs$data$pangenomes[[pangenome]]$files$metadata, genus = confs$genus) %>% 
  dplyr::mutate(
    virulence_pcr = dplyr::case_when(
      virulence_pcr == "positive" ~ "+",
      virulence_pcr == "negative" ~ "-",
      TRUE ~ virulence_pcr
    )
  )

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

probeInfo <- suppressMessages(readr::read_tsv(confs$data$other$files$probe_info)) %>% 
  dplyr::mutate(
    type = forcats::fct_relevel(.f = type, "forward", "probe", "reverse")
  )

################################################################################
## extract the ANI tree information for the node order
rawTree <- ape::read.tree(file = confs$analysis$phylogeny$ani_upgma$files$tree_rooted)

## set negative length edges => 0
rawTree$edge.length[rawTree$edge.length < 0] <- 0


## P. brasiliense species clade
nodeOfInterest <- dplyr::filter(sampleInfo, SpeciesName == "P. brasiliense") %>% 
  dplyr::pull(Genome)

clade <- ape::getMRCA(phy = rawTree, tip = nodeOfInterest)
subTree <- ape::extract.clade(phy = rawTree, node = clade)

sampleInfo <- dplyr::left_join(
  x = tibble::tibble(Genome = nodeOfInterest),
  y = sampleInfo, by = "Genome"
)

## tree for plotting 
## add data to tree
treeTbl <- as_tibble(subTree) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata()

pt_treeMain <- ggtree::ggtree(
  tr = treeTbl
)

leafOrder <- dplyr::arrange(.data = pt_treeMain$data, y) %>% 
  dplyr::filter(isTip) %>% 
  dplyr::pull(label)

# leafOrder <- rev(ggtree::get_taxa_name(pt_treeMain))

sampleInfo %<>% 
  dplyr::mutate(
    Genome = forcats::fct_relevel(Genome, !!!leafOrder)
  )

################################################################################
## process assay blast result
blastRes <- suppressMessages(
  readr::read_tsv(
    file = confs$analysis$insilico_assay$files$blastn,
    col_names = c(
      "qseqid", "qstart", "qend", "qlen", "sseqid", "sstart", "send", "sstrand", "slen",
      "pident", "length", "mismatch", "qcovs", "gapopen", "evalue", "bitscore"
    )
  )
) %>% 
  dplyr::mutate(
    sstrand = dplyr::case_when(
      sstrand == "plus" ~ "+",
      sstrand == "minus" ~ "-",
    )
  )


blastGr <- GenomicRanges::makeGRangesFromDataFrame(
  df = dplyr::mutate(
    .data = blastRes,
    chr = sseqid,
    start = dplyr::if_else(sstrand == "+", true = sstart, false = send),
    end = dplyr::if_else(sstrand == "+", true = send, false = sstart),
    strand = sstrand
  ) %>% 
    dplyr::select(chr, start, end, strand, everything()),
  keep.extra.columns = TRUE
)

## extend by 100bp and merge to generate the broad GRanges covering all PCR regions (F, P, R)
reducedGr <- GenomicRanges::reduce(
  x = blastGr, min.gapwidth = 100, ignore.strand = TRUE
)

mcols(reducedGr)$regionId <- paste("region_", 1:length(reducedGr), sep = "")

ovlpHits <- GenomicRanges::findOverlaps(
  query = blastGr, subject = reducedGr
)

if(any(duplicated(ovlpHits@from))){
  stop("Duplicated hits found in GRanges findOverlaps()")
}

mcols(blastGr)$regionId <- NA
mcols(blastGr)$regionId[ovlpHits@from] <- mcols(reducedGr)$regionId[ovlpHits@to]

assayBlastDf <- dplyr::left_join(
  x =  tibble::as_tibble(mcols(blastGr)), y = genomeChrs, by = c("sseqid" = "chr")
) %>% 
  # dplyr::left_join(
  #   y = sampleInfo, by = "Genome"
  # ) %>% 
  dplyr::left_join(
    y = probeInfo, by = c("qseqid" = "probeId")
  ) %>% 
  dplyr::group_by(Genome, assay, sseqid, regionId) %>% 
  dplyr::mutate(
    assayRegionHits = n(),
    assayGenomeGrp = cur_group_id(),
    assayMismatches = sum(mismatch)
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Genome, assay, sseqid, sstart, type) %>% 
  dplyr::select(Genome,
                # sampleId, SpeciesName, strain,
                assay, type,
                assayRegionHits, regionId, assayGenomeGrp, everything())


readr::write_tsv(x = assayBlastDf, file = paste(outPrefix, ".data.tab", sep = ""))

## select only best matching assay
bestAssaysInGenomes <- dplyr::filter(assayBlastDf, assay != "P_wasabiae_WPP163") %>% 
  dplyr::distinct(
    assayGenomeGrp, Genome, assay, sseqid, regionId, assayRegionHits, assayMismatches
  ) %>% 
  dplyr::group_by(Genome, assay) %>% 
  dplyr::arrange(desc(assayRegionHits), assayMismatches, .by_group = TRUE) %>% 
  dplyr::slice(1L) %>% 
  dplyr::ungroup()

bestAssaysInGenomes <- dplyr::select(bestAssaysInGenomes, assayGenomeGrp) %>% 
  dplyr::left_join(y = assayBlastDf, by = "assayGenomeGrp") %>% 
  dplyr::filter(assayRegionHits == 3) %>% 
  dplyr::mutate(
    primerCov = if_else(condition = qcovs == 100, true = "complete", false = "incomplete"),
    # mismatch = if_else(condition = mismatch == 0, true = "no", false = "yes"),
    primerCov = forcats::as_factor(primerCov)
  )

## finally select assays only for the Genomes in sampleInfo
genomeAssayDetails <- dplyr::select(sampleInfo, Genome) %>% 
  dplyr::left_join(y = bestAssaysInGenomes, by = "Genome") %>% 
  dplyr::filter(!is.na(assay)) %>% 
  dplyr::mutate(
    Genome = forcats::fct_relevel(Genome, !!!leafOrder),
    assay = forcats::fct_relevel(
      assay, "Erwinia", "P_br_1360500", "LZI", "TIR"
    )
  )

genomeAssaySummary <- dplyr::group_by(
  bestAssaysInGenomes, assayGenomeGrp, Genome, assay, assayMismatches
) %>% 
  dplyr::summarise(
    assayLen = sum(qlen), .groups = "drop"
  ) %>% 
  dplyr::mutate(
    assayPident = 100*(assayLen - assayMismatches)/assayLen,
    ####************************************************************************
    ####*This calculation is not correct
    ####*Need to use the original assay length from FASTA file
    ####*Because, if `dplyr::filter(assayRegionHits == 3)` is changed to 
    ####*include assays where 2 of F,P,R mapped, assayLen calculation will be
    ####*wrong as it will exclude the missing probe length.
    ####*Solution: update the probe metadata file to include the length too.
    ####************************************************************************
    
    assayResult = if_else(assayMismatches == 0, "pass", "fail")
  )

################################################################################

pt_tree <- pt_treeMain

## assay blastn visualization
(
  pt_pcrBlast <- ggplot2::ggplot(
    data = genomeAssayDetails,
    mapping = aes(x = type, y = Genome)
  ) +
    geom_tile(
      mapping = aes(fill = pident, width = qcovs/100, height = qcovs/100),
      color = "black", linewidth = 0.5) +
    labs(
      title = "blastn results for forward (F), probe (P) and reverse (R) PCR primers"
    ) +
    facet_grid(cols = vars(assay), scales = "free_y") +
    viridis::scale_fill_viridis(name = "% identity", option = "viridis", direction = -1) +
    scale_x_discrete(
      labels = c("forward" = "F", "probe" = "P", "reverse" = "R")
    ) +
    guides(
      size = guide_legend(override.aes = list(fill = "white"))
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14, face = "bold"),
      plot.title = element_text(hjust = 0.5, vjust = 0),
      strip.text.x = element_text(size = 12, face = "bold") 
    )
)

##################################################################################
## generate annotation plots
pt_theme <- theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, face = "bold", vjust = 0.5),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank()
  )

pt_source <- dplyr::select(sampleInfo, Genome, source) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "source", y = Genome, color = source)
  ) +
  geom_point(shape = 15, size = 1) +
  scale_color_viridis_d() +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )

pt_region <- dplyr::select(sampleInfo, Genome, geo_loc_country) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "country", y = Genome, color = geo_loc_country)
  ) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("Netherlands" = "blue")
  ) +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )

pt_vir <- dplyr::select(sampleInfo, Genome, virulence) %>% 
  dplyr::mutate(
    virulence = forcats::fct_relevel(virulence, "virulent", "avirulent")
  ) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "virulence", y = Genome, shape = virulence, color = virulence)
  ) +
  geom_point(size = 1) +
  scale_shape_manual(
    name = "Virulence",
    values = c("virulent" = 17, "avirulent" = 19)
  ) +
  scale_color_manual(
    name = "Virulence",
    values = c("virulent" = "red", "avirulent" = "green")
  ) +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )


pt_pcr <- dplyr::select(sampleInfo, Genome, virulence_pcr) %>% 
  ggplot2::ggplot(
    mapping = aes(x = "PCR", y = Genome, label = virulence_pcr, color = virulence_pcr)
  ) +
  geom_text(size = 4, fontface = "bold") +
  scale_color_manual(
    values = c("+" = "red", "-" = "green"), breaks = NULL
  ) +
  pt_theme +
  theme(
    axis.ticks = element_blank()
  )

pt_sampleYear <- dplyr::select(sampleInfo, Genome, collection_year, source) %>%
  ggplot2::ggplot(
    mapping = aes(x = collection_year, y = Genome, color = source)
  ) +
  geom_point(size = 1) +
  scale_color_manual(
    values = c("NCBI" = "black", "NAK" = "#468e30", "NVWA" = "black"),
    breaks = NULL
  ) +
  scale_x_continuous(expand = expansion(add = 0.5)) +
  coord_cartesian(xlim = c(2000, 2022)) +
  labs(title = "Sampling\nyear") +
  pt_theme +
  theme(
    plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0)
  )


## arrange plots one by one
pt_all <- pt_source %>% aplot::insert_left(pt_tree, width = 12) %>%
  aplot::insert_right(pt_sampleYear, width = 8) %>%
  aplot::insert_right(pt_region, width = 1) %>%
  aplot::insert_right(pt_vir, width = 1) %>%
  aplot::insert_right(pt_pcr, width = 0.5) %>%
  aplot::insert_right(pt_pcrBlast, width = 15)


# png(filename = paste(outPrefix, ".png", sep = ""), width = 3000, height = 1500, res = 200)
pdf(file = paste(outPrefix, ".pdf", sep = ""), width = 15, height = 12)
pt_all
dev.off()


pt_treeBlast <- pt_vir %>% aplot::insert_right(pt_tree, width = 12) %>%
  aplot::insert_right(pt_pcrBlast, width = 8)


# png(filename = paste(outPrefix, ".treeBlast.png", sep = ""), width = 2500, height = 1500, res = 150)
pdf(file = paste(outPrefix, ".treeBlast.pdf", sep = ""), width = 12, height = 12)
pt_treeBlast
dev.off()

