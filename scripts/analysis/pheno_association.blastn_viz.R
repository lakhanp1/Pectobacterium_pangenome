#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(configr))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(treeio))
suppressPackageStartupMessages(library(ggtree))
suppressPackageStartupMessages(library(openxlsx))

## parse blast performed against pangenome for phenotype specific sqeuences
## visualize the results
## write the gene ids that show no blast hits outside the genome of interest

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/association_analysis.R")
source("scripts/utils/phylogeny_functions.R")
################################################################################
set.seed(124)

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
treeMethod <- "kmer_nj" # ani_upgma, kmer_nj
outGroup <- confs$analysis$phylogeny$outgroup

phenotype <- "assay_FN"
file_blastOut <- file.path(
  confs$analysis$association$path, phenotype, "pheno_specific.pangenome.blastn.out"
)

outDir <- file.path(confs$analysis$association$path, phenotype)
outPrefix <- file.path(outDir, phenotype)

file_associatedSeqInfo <- paste(outPrefix, ".pheno_specific.seq_info.txt", sep = "")

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, Genome
)

rawTree <- ape::read.tree(file = confs$analysis$phylogeny[[treeMethod]]$files$tree)

genomeChrs <- suppressMessages(
  readr::read_tsv(
    file = panConf$files$genome_chrs,
    col_names = c("Genome", "chr"), col_types = "c"
  )
)

## genomes that are used for phenotype association analysis
associatedGenomes <- get_phenotype_association_genomes(
  phenotype = phenotype,
  confFile = panConf$analysis_confs$files$phenotype_association
)

## sequence info: use for chr wise separation and ordering
hgSeqInfo <- suppressMessages(readr::read_tsv(file_associatedSeqInfo)) %>%
  dplyr::mutate(
    dplyr::across(.cols = c(Genome, homology_group_id), .fns = ~ as.character(.x))
  ) %>%
  dplyr::select(
    qGenome = Genome, qseqid = mRNA_identifier, qChr = chr, qChrStart = start,
    qChrEnd = end, qHomologyGrp = homology_group_id
  ) %>%
  dplyr::arrange(qGenome, qChr, qChrStart)

## import blast result
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
    totalqcov = as.numeric(sprintf(fmt = "%.2f", 100 * (qend - qstart + 1) / qlen))
  ) %>%
  dplyr::left_join(y = hgSeqInfo, by = "qseqid") %>%
  dplyr::left_join(
    y = dplyr::rename(genomeChrs, sGenome = Genome), by = c("sseqid" = "chr")
  ) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, Genome, sSpecies = SpeciesName, sSample = sampleName),
    by = c("sGenome" = "Genome")
  ) %>%
  dplyr::left_join(
    ## add tag to blast hit if it is from the phenotype specific genome
    y = tibble::enframe(associatedGenomes, name = "sGenomeType", value = "sGenome") %>%
      tidyr::unnest(cols = sGenome),
    by = "sGenome"
  ) %>%
  tidyr::replace_na(list(sGenomeType = "other")) %>%
  dplyr::relocate(sGenome, sGenomeType, sSpecies, sSample, .before = sseqid) %>%
  dplyr::relocate(qGenome, qChr, qChrStart, qChrEnd, qHomologyGrp, .after = qseqid) %>%
  dplyr::relocate(totalqcov, .before = qcovs)

################################################################################
## prepare data for visualization

## root tree
rootedTr <- ape::root(
  phy = rawTree, outgroup = sampleInfoList[[outGroup]]$Genome, edgelabel = TRUE
) %>%
  ape::ladderize()

## add data to tree
treeTbl <- as_tibble(rootedTr) %>%
  dplyr::full_join(y = sampleInfo, by = c("label" = "Genome")) %>%
  treeio::as.treedata()

pt_tree <- ggtree::ggtree(tr = treeTbl) +
  labs(title = confs$analysis$phylogeny[[treeMethod]]$name) +
  ggtree::geom_tiplab(
    mapping = aes(label = nodeLabs),
    size = 3, align = TRUE, linesize = 0.5
  )

leafOrder <- ggtree::get_taxa_name(pt_tree)

## filter blastn results
significantHits <- dplyr::filter(
  blastRes, totalqcov >= 80, pident >= 70
)

bestHits <- dplyr::group_by(significantHits, qseqid, sGenome) %>%
  dplyr::arrange(desc(totalqcov), desc(pident), .by_group = TRUE) %>%
  dplyr::slice(1L) %>%
  dplyr::ungroup() %>%
  dplyr::select(qseqid, qHomologyGrp, qGenome, qChr, qChrStart, sGenome, pident, sGenomeType) %>%
  dplyr::arrange(qGenome, qChr, qChrStart)

## ensure that the genomes without any blast hits are added for proper alignment with tree
bestHits %<>%
  dplyr::full_join(y = tibble::tibble(sGenome = leafOrder), by = "sGenome") %>%
  tidyr::replace_na(list(qseqid = "genome-without-match")) %>%
  dplyr::mutate(
    qseqid = forcats::as_factor(qseqid)
  )

## genes that matches only in compare genome category and not in any other
specificAssays <- setdiff(
  x = bestHits$qseqid,
  y = dplyr::filter(bestHits, sGenomeType != "compare") %>%
    dplyr::distinct(as.character(qseqid)) %>% tibble::deframe() %>%
    append("genome-without-match")
)

tibble::tibble(seqId = specificAssays) %>%
  readr::write_tsv(file = paste(outPrefix, ".specific_assay_seq_id.txt", sep = ""))

################################################################################
## visualize result along with evolutionary tree

## mark outgroup
pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outGroup, column = "sampleName")

## mark species of interest
pt_tree3 <- pt_tree2 +
  ggtree::geom_tippoint(
    mapping = aes(subset = c(label %in% c(associatedGenomes$compare))),
    color = "#088734"
  ) +
  ggtree::geom_tippoint(
    mapping = aes(subset = c(label %in% c(associatedGenomes$against))),
    color = "red"
  ) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  ggnewscale::new_scale_color() +
  ## virulence phenotype
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "virulence", color = virulence),
    geom = "geom_point", shape = 17, size = 2,
    pwidth = 0.01, offset = 0.2
  ) +
  scale_color_manual(
    values = c("virulent" = "red", "avirulent" = "green"),
    na.value = alpha("white", 0)
  ) +
  ggnewscale::new_scale_color() +
  ## virulence PCR result
  ggtreeExtra::geom_fruit(
    mapping = aes(y = id, x = "vir_pcr", color = virulence_pcr),
    geom = "geom_point",
    pwidth = 0.01, offset = 0.01
  ) +
  scale_color_manual(
    values = c("positive" = "red", "negative" = "green"),
    na.value = alpha("white", 0)
  )


## blastn visualization
pt_blast <- ggplot2::ggplot(
  data = bestHits,
  mapping = aes(x = qseqid, y = sGenome)
) +
  geom_tile(
    mapping = aes(fill = pident),
    color = "black", linewidth = 0.01
  ) +
  geom_vline(
    data = tibble(qseqid = specificAssays),
    mapping = aes(xintercept = qseqid), color = "green"
  ) +
  labs(
    title = "blastn results %identity"
  ) +
  scale_fill_viridis_c(name = "% identity", option = "plasma") +
  scale_x_discrete(expand = expansion(add = 0)) +
  facet_grid(cols = vars(qChr), scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(
    panel.grid.major.x = element_blank(),
    axis.title = element_blank(),
    # axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5),
    plot.title = element_text(hjust = 0.5, vjust = 0),
    # strip.text.x = element_text(size = 10, face = "bold"),
    strip.text = element_blank()
  )

## arrange plots one by one
pt_all <- pt_blast %>% aplot::insert_left(pt_tree3, width = 1)

ggsave(
  plot = pt_all, width = 14, height = 20, scale = 2,
  filename = paste(outPrefix, ".specific_seq_blast.pdf", sep = "")
)

################################################################################
## write blast output
wb <- openxlsx::createWorkbook()

currentSheet <- "blast_result"
openxlsx::addWorksheet(wb, sheetName = currentSheet)
openxlsx::writeDataTable(
  wb = wb, sheet = currentSheet, x = blastRes, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)

openxlsx::freezePane(wb = wb, sheet = currentSheet, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
openxlsx::saveWorkbook(
  wb = wb,
  file = file.path(outDir, "blast_result.xlsx"), overwrite = TRUE
)

################################################################################
