suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(GenomicRanges))

## parse PCR primer blastn results and generate summary figure

rm(list = ls())


################################################################################

analysisName <- "bioassay_probe_match"

outDir <- here::here("analysis", "04_pangenome.50g", "primer_blastn")
outPrefix <- file.path(outDir, analysisName)


file_blastOut <- here::here("analysis", "04_pangenome.50g", "primer_blastn", "probe_blasnt.out")
file_phenotypes <- here::here("analysis", "04_pangenome.50g", "pectobacterium.50g.phenotypes.txt")
file_genomeChr <- here::here("analysis", "04_pangenome.50g", "primer_blastn", "genome_chr.tab")
file_probeInfo <- here::here("analysis", "04_pangenome.50g", "primer_blastn", "probe_info.txt")

pangenome <- here::here("analysis", "04_pangenome.50g", "pectobacterium.50g.DB")
file_aniTree <- file.path(pangenome, "ANI", "fastANI", "ANI_RENAMED.newick")

################################################################################

genomeChrs <- suppressMessages(
  readr::read_tsv(file = file_genomeChr, col_names = c("genomeId", "chr"))
)
phenotypes <- suppressMessages(readr::read_tsv(file_phenotypes)) %>% 
  tidyr::unite(genomeId, Genome, id, sep = "_", remove = F)

probeInfo <- suppressMessages(readr::read_tsv(file_probeInfo)) %>% 
  dplyr::mutate(
    type = forcats::fct_relevel(.f = type, "forward", "probe", "reverse")
  )

################################################################################
## extract the ANI tree information for the node order
aniTreeTibble <- treeio::read.newick(file = file_aniTree) %>% 
  # treeio::as.treedata() %>% 
  as_tibble()

aniTree <- dplyr::full_join(
  x = aniTreeTibble, y = phenotypes, by = c("label" = "genomeId")
) %>% 
  treeio::as.treedata()


pt_tree <- ggtree::ggtree(
  tr = aniTree,
  # mapping = aes(color = species)
  # branch.length = "none"
) +
  scale_x_continuous(
    expand = expansion(add = c(0.01, 0.06))
  ) +
  ggtree::geom_tiplab(
    mapping = aes(color = species, label = id),
    align = TRUE
  ) +
  ggtree::geom_treescale(x = 0, y = length(aniTree@phylo$tip.label)-2) 

leafOrder <- dplyr::arrange(.data = pt_tree$data, y) %>% 
  dplyr::filter(isTip) %>% 
  dplyr::pull(label)

phenotypes <- dplyr::mutate(
  .data = phenotypes,
  genomeId = forcats::fct_relevel(genomeId, !!!leafOrder),
  pcr = dplyr::case_when(
    pcr == "positive" ~ "+",
    pcr == "negative" ~ "-",
    TRUE ~ pcr
  )
)

################################################################################

blastRes <- suppressMessages(
  readr::read_tsv(
    file = file_blastOut,
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

mergedDf <- dplyr::left_join(
  x =  tibble::as_tibble(mcols(blastGr)), y = genomeChrs, by = c("sseqid" = "chr")
) %>% 
  dplyr::left_join(
    y = phenotypes, by = "genomeId"
  ) %>% 
  dplyr::left_join(
    y = probeInfo, by = c("qseqid" = "probeId")
  ) %>% 
  dplyr::group_by(Genome, assay, sseqid, regionId) %>% 
  dplyr::mutate(
    n = n()
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::arrange(Genome, assay, sseqid, sstart, type) %>% 
  dplyr::mutate(
    genomeId = forcats::fct_relevel(genomeId, !!!leafOrder)
  ) %>% 
  dplyr::select(Genome, genomeId, strain, species, assay, type, n, everything())


readr::write_tsv(x = mergedDf, file = paste(outPrefix, ".data.tab", sep = ""))



pcrPlotDf <- dplyr::filter(mergedDf, assay != "P_wasabiae_WPP163") %>% 
  dplyr::mutate(
    primerCov = if_else(condition = qcovs == 100, true = "complete", false = "incomplete"),
    mismatch = if_else(condition = mismatch == 0, true = "no", false = "yes")
  ) %>% 
  dplyr::filter(n != 1)

table(pcrPlotDf$mismatch, pcrPlotDf$n)

pt_pcrBlast <- ggplot2::ggplot(
  data = pcrPlotDf,
  mapping = aes(x = type, y = genomeId)
) +
  geom_tile(
    mapping = aes(fill = pident, size = primerCov,
               width = qcovs/100, height = qcovs/100),
    color = "black") +
  geom_point(mapping = aes(alpha = mismatch), size = 2) +
  labs(
    title = "blastn results for forward (F), probe (P) and reverse (R) PCR primers"
  ) +
  scale_fill_viridis_b(name = "% identity", option = "plasma") +
  scale_alpha_manual(
    name = NULL,
    values = c("no" = 0, "yes" = 1),
    breaks = c("yes"),
    labels = c("mismatch")
  ) +
  scale_size_manual(
    name = "Primer coverage",
    values = c("complete" = 1.5, "incomplete" = 0.5)
  ) +
  scale_x_discrete(
    labels = c("forward" = "F", "probe" = "P", "reverse" = "R")
  ) +
  guides(
    size = guide_legend(override.aes = list(fill = "white"))
  ) +
  facet_grid(cols = vars(assay, n), scales = "free_y") +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    plot.title = element_text(hjust = 0.5, vjust = 0)
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
  scale_fill_viridis_c(
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

# pt_region <- dplyr::select(phenotypes, genomeId, region) %>% 
#   ggplot2::ggplot(
#     mapping = aes(x = "Region", y = genomeId, label = region)
#   ) +
#   geom_text(size = 4) +
#   labs(title = "Region") +
#   pt_theme +
#   theme(
#     plot.title = element_text(angle = 0, hjust = 0.5, vjust = 0),
#     axis.text.x = element_blank(),
#     axis.ticks = element_blank()
#   )

## arrange plots one by one
pt_all <- pt_vir %>% aplot::insert_left(pt_tree, width = 15) %>%
  aplot::insert_right(pt_pcr, width = 0.5) %>%
  aplot::insert_right(pt_symptoms, width = 1) %>%
  aplot::insert_right(pt_sampleYear, width = 2) %>%
  aplot::insert_right(pt_tissue, width = 1.5) %>%
  aplot::insert_right(pt_pcrBlast, width = 12)


png(filename = paste(outPrefix, ".23.png", sep = ""), width = 3000, height = 1500, res = 150)
pt_all
dev.off()







