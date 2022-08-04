suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(GenomicRanges))

## parse PCR primer blastn results and generate summary figure

rm(list = ls())


################################################################################

analysisName <- "LZI_TIR_gene_blastn"

outDir <- here::here("analysis", "04_pangenome.50g", "primer_blastn")
outPrefix <- file.path(outDir, analysisName)


file_blastOut <- here::here("analysis", "04_pangenome.50g", "primer_blastn", "LZI_TIR_blastn.out")
file_phenotypes <- here::here("analysis", "04_pangenome.50g", "pectobacterium.50g.phenotypes.txt")
file_genomeChr <- here::here("analysis", "04_pangenome.50g", "primer_blastn", "genome_chr.tab")

pangenome <- here::here("analysis", "04_pangenome.50g", "pectobacterium.50g.DB")
file_aniTree <- file.path(pangenome, "ANI", "fastANI", "ANI_RENAMED.newick")

################################################################################
genomeChrs <- suppressMessages(
  readr::read_tsv(file = file_genomeChr, col_names = c("genomeId", "chr"))
)
phenotypes <- suppressMessages(readr::read_tsv(file_phenotypes)) %>% 
  tidyr::unite(genomeId, Genome, id, sep = "_", remove = F)

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


mergedDf <- dplyr::left_join(
  x =  blastRes, y = genomeChrs, by = c("sseqid" = "chr")
) %>% 
  dplyr::left_join(
    y = phenotypes, by = "genomeId"
  ) %>% 
  dplyr::group_by(Genome, sseqid) %>%
  dplyr::mutate(
    n = n()
  ) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(qseqid, Genome, sseqid, sstart) %>% 
  dplyr::select(Genome, genomeId, strain, species, n, everything())


readr::write_tsv(x = mergedDf, file = paste(outPrefix, ".data.tab", sep = ""))







