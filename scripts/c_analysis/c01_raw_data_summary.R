suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(here))
suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(XML))
suppressPackageStartupMessages(library(openxlsx))


# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
# Script to summarize raw data and generate reports and tables

rm(list = ls())

#####################################################################
analysisName <- "raw_data_summary"
outDir <- here::here("analysis", "02_raw_data_summary")

file_assembly <- here::here("data/reference_data", "assembly_docsum.xml")
file_biosample <- here::here("data/reference_data", "biosample_docsum.xml")
file_inhouseMmetadata <- here::here("data/reference_data", "inhouse_samples_metadata.txt")
file_buscoMqc <- here::here("analysis/01_multiqc", "busco_multiqc_data/multiqc_busco.txt")
file_quastMqc <- here::here("analysis/01_multiqc", "quast_multiqc_data/multiqc_quast.txt")

#####################################################################
## inhouse genomes and QC data
inhouseGenomes <- suppressMessages(readr::read_tsv(file = file_inhouseMmetadata)) %>% 
  dplyr::mutate(
    source = "inhouse",
    SubmissionDate = "2023"
  )

buscoMqc <- suppressMessages(readr::read_tsv(file = file_buscoMqc)) %>% 
  dplyr::mutate(
    sampleId = stringr::str_replace(
      string = Sample, pattern = "short_summary.specific.enterobacterales_odb10.",
      replacement = ""
    )
  ) %>% 
  dplyr::select(-Sample, -lineage_dataset)

quastCols <- c(
  contigs_gt_0kb = "# contigs (>= 0 bp)",
  contigs_gt_1kb = "# contigs (>= 1000 bp)",
  contigs_gt_5kb = "# contigs (>= 5000 bp)",
  contigs_gt_10kb = "# contigs (>= 10000 bp)",
  contigs_gt_25kb = "# contigs (>= 25000 bp)",
  contigs_gt_50kb = "# contigs (>= 50000 bp)",
  length_gt_0kb = "Total length (>= 0 bp)",
  length_gt_1kb = "Total length (>= 1000 bp)",
  length_gt_5kb = "Total length (>= 5000 bp)",
  length_gt_10kb = "Total length (>= 10000 bp)",
  length_gt_25kb = "Total length (>= 25000 bp)",
  length_gt_50kb = "Total length (>= 50000 bp)",
  n_contigs = "# contigs",
  largest_contig = "Largest contig",
  length = "Total length",
  GC_per = "GC (%)",
  N50 = "N50",
  N90 = "N90",
  auN = "auN",
  L50 = "L50",
  L90 = "L90",
  N_per_100kb = "# N's per 100 kbp"
)

quastMqc <- suppressMessages(readr::read_tsv(file = file_quastMqc)) %>% 
  dplyr::rename(
    sampleId = Sample,
    !!!quastCols
  )

#####################################################################
## extract assembly information
asmDoc <- XML::xmlParse(file = file_assembly)
## asmDoc <- XML::xmlTreeParse(file = file_assembly)
asmNodes <- XML::getNodeSet(doc = asmDoc, path = "//DocumentSummary")
asmDf <- XML::xmlToDataFrame(nodes = asmNodes)

# xmlList <- XML::xmlToList(node = asmNodes[[1]])
# intersect(unlist(xmlList$PropertyList), c("latest_genbank", "latest_refseq"))

asmMetadata <- purrr::map_dfr(
  .x = asmNodes,
  .f = function(x){
    xmlList <- XML::xmlToList(node = x)
    return(c(
      Id = xmlList$Id,
      taxonomy_check_status = xmlList$Meta$`taxonomy-check-status`,
      representative_status = xmlList$Meta$`representative-status`,
      RS_BioprojectAccn = xmlList$RS_BioProjects$Bioproj$BioprojectAccn,
      GB_BioprojectAccn = xmlList$GB_BioProjects$Bioproj$BioprojectAccn
    ))
  }
)

#####################################################################
## extract BioSample information
bsDoc <- XML::xmlParse(file = file_biosample)
bsNodes <- XML::getNodeSet(doc = bsDoc, path = "//DocumentSummary")

colInfo <- NULL

# aNode <- bsNodes[[1]]
# aNodeList <- XML::xmlToList(node = aNode)

bsMetadata <- purrr::map_dfr(
  .x = bsNodes,
  .f = function(aNode){
    aNodeDoc <- XML::xmlParse(XML::toString.XMLNode(aNode))
    aNodeList <- XML::xmlToList(node = aNode)
    
    ## parse all the metadata from Attribute nodes
    attrDf <- purrr::map_dfr(
      .x = XML::getNodeSet(doc = aNodeDoc, path = "//Attribute"),
      .f = function(nd){
        return(c(
          XML::xmlAttrs(node = nd),
          value = XML::xmlValue(x = nd)
        ))
      }
    )
    
    ## cleanup attribute names which will be used as column names
    attrDf <- dplyr::mutate(
      .data = attrDf,
      attribute_name = tolower(
        stringr::str_replace_all(
          string = attribute_name, pattern = "(\\W+)", replacement = "_"
        )
      ),
      attribute_name = stringr::str_replace(
        string = attribute_name, pattern = "_$", replacement = ""
      ),
      harmonized_name = dplyr::coalesce(harmonized_name, attribute_name)
    )
    
    ## cleanup the attribute values
    attrDf <- dplyr::mutate(
      .data = attrDf,
      value = stringr::str_replace_all(
        string = value,
        pattern = regex(
          pattern = "(not.applicable|missing|unknown|not.provided)",
          ignore_case = TRUE
        ),
        replacement = NA_character_
      )
    )
    
    ## add the attribute name and description to global DF
    colInfo <<- dplyr::bind_rows(
      colInfo,
      dplyr::filter(attrDf, !is.na(value)) %>% 
        dplyr::select(attribute_name, harmonized_name, display_name)
    )
    
    return(c(
      BioSampleAccn = aNodeList$Accession,
      tibble::deframe(x = dplyr::select(attrDf, harmonized_name, value))
    ))
  }
)

colStats <- dplyr::count(colInfo, harmonized_name, display_name, sort = TRUE)

bsMetadata <- dplyr::select(
  bsMetadata, BioSampleAccn, colStats$harmonized_name
)


#####################################################################
## merge all data
ncbiData <- dplyr::select(
  .data = asmDf,
  AssemblyAccession, AssemblyName, SpeciesName, SpeciesTaxid, Organism, 
  AssemblyStatus, AssemblyType, BioSampleAccn, SubmissionDate, SubmitterOrganization,
  FtpPath_RefSeq, FtpPath_GenBank, assembly_id = Id
) %>% 
  dplyr::left_join(y = asmMetadata, by = c("assembly_id" = "Id")) %>%
  dplyr::left_join(y = bsMetadata, by = "BioSampleAccn") %>% 
  dplyr::mutate(
    AssemblyName = stringr::str_replace_all(
      string = AssemblyName, pattern = "\\s+", replacement = "_"
    ),
    source = "NCBI"
  ) %>% 
  tidyr::unite(col = sampleId, AssemblyAccession, AssemblyName, sep = "_", remove = F)


genomeMetadata <- dplyr::bind_rows(
  ncbiData, inhouseGenomes
) %>% 
  dplyr::left_join(y = buscoMqc, by = "sampleId") %>% 
  dplyr::left_join(y = quastMqc, by = "sampleId")

genomeMetadata <- dplyr::mutate(
  .data = genomeMetadata,  
  SubmissionDate = lubridate::ymd_hm(SubmissionDate),
  submission_y = lubridate::year(SubmissionDate),
  submission_m = lubridate::month(SubmissionDate)
)


readr::write_tsv(
  x = genomeMetadata,
  file = here::here("data/reference_data", "sample_metadata.tsv")
)

wb <- openxlsx::createWorkbook()
openxlsx::addWorksheet(wb, sheetName = "metadata")
openxlsx::writeDataTable(
  wb = wb, sheet = 1, x = genomeMetadata, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 1, firstActiveRow = 2, firstActiveCol = 2)
openxlsx::addWorksheet(wb, sheetName = "column_info")
openxlsx::writeDataTable(
  wb = wb, sheet = 2, x = colStats, withFilter = TRUE,
  keepNA = TRUE, na.string = "NA"
)
openxlsx::freezePane(wb = wb, sheet = 2, firstActiveRow = 2, firstActiveCol = 2)

# openxlsx::openXL(wb)
saveWorkbook(
  wb = wb, file = here::here("data/reference_data", "sample_metadata.xlsx"), overwrite = TRUE
)

#####################################################################
pt_theme <- theme_bw() +
  theme(
    panel.grid = element_blank(),
    axis.text = element_text(face = "bold", size = 16, color = "black"),
    plot.title = element_text(face = "bold", size = 18, color = "black"),
    axis.title = element_text(face = "bold", size = 16, color = "black")
  )

pt_speciesCount <- ggplot2::ggplot(
  data = genomeMetadata,
  mapping = aes(y = forcats::fct_infreq(SpeciesName))
) +
  geom_bar(fill = "black") +
  labs(
    x = "#genomes",
    y = NULL,
    title = "# of genomes"
  ) +
  scale_x_continuous(expand = expansion(add = c(0, 2))) +
  pt_theme +
  theme(
    axis.text.y = element_text(face = "bold.italic", size = 14, color = "black"),
    axis.text.x = element_text(face = "bold", size = 16, color = "black"),
    plot.title = element_text(face = "bold", size = 18, color = "black"),
    axis.title = element_text(face = "bold", size = 16, color = "black")
  )

ggsave(
  filename = file.path(outDir, "species_count_bar.pdf"),
  plot = pt_speciesCount, width = 10, height = 8
)


#####################################################################

summaryStatsDf <- pivot_longer(
  data = dplyr::select(
    genomeMetadata, sampleId, AssemblyStatus, taxonomy_check_status, source
  ),
  cols = -sampleId,
  names_to = "variable"
) %>% 
  tidyr::replace_na(
    replace = list(value = "NA")
  ) %>% 
  dplyr::count(variable, value, name = "count") %>% 
  tidyr::unite(col = "label", value, count, sep = ": ", remove = F)

pt_otherSummary <- ggplot2::ggplot(
  data = summaryStatsDf,
  mapping = aes(x = variable, y = count, fill = value)
) +
  geom_bar(
    stat = "identity", color = "grey", alpha = 0, size = 1
  ) +
  geom_label(
    mapping = aes(label = label),
    position = position_stack(vjust = 0.5), color = "black", alpha = 0, label.size = 1
  ) +
  pt_theme +
  theme(
    legend.position = "none",
    axis.text.x = element_text(face = "bold", size = 16, color = "black"),
    plot.title = element_text(face = "bold", size = 18, color = "black"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()
  )


ggsave(
  filename = file.path(outDir, "other_summary_bar.pdf"),
  plot = pt_otherSummary, width = 7, height = 6
)

#####################################################################

## n_contigs histogram
ptHist_contigs <- ggplot2::ggplot(
  data = genomeMetadata,
  mapping = aes(x = n_contigs)
) +
  geom_histogram(bins = 100, color = "black", fill = "black") + 
  labs(
    title = "Histogram of contig numbers",
    x = "# contigs", y = "Count"
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "N_contigs.hist.pdf"),
  plot = ptHist_contigs, width = 8, height = 8
)


## n_contigs histogram excluding outliers
ptHist_contigs2 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = n_contigs)
) +
  geom_histogram(bins = 100, color = "black", fill = "black") + 
  labs(
    title = "Histogram of contig numbers (excluding 2 outliers)",
    x = "# contigs", y = "Count"
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "N_contigs2.hist.pdf"),
  plot = ptHist_contigs2, width = 8, height = 8
)


## N50 histogram
ptHist_n50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = N50)
) +
  geom_histogram(bins = 100, color = "black", fill = "black") + 
  scale_x_continuous(
    labels = label_comma(scale_cut = c(bp = 0, kb = 1000, mb = 1000000))
  ) +
  labs(
    title = "Histogram of N50 (excluding 2 outliers)",
    x = "N50", y = "Count"
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "N50.hist.pdf"),
  plot = ptHist_n50, width = 8, height = 8
)


ptHist_l50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = L50)
) +
  geom_histogram(bins = 100, color = "black", fill = "black") + 
  labs(
    title = "Histogram of L50 (excluding 2 outliers)",
    x = "L50", y = "Count"
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "L50.hist.pdf"),
  plot = ptHist_l50, width = 8, height = 8
)


## n_contigs vs N50
ptScat_contig_vs_n50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = N50, y = n_contigs)
) +
  geom_point() +
  labs(title = "#contigs vs N50") +
  scale_x_continuous(
    labels = label_comma(scale_cut = c(bp = 0, kb = 1000, mb = 1000000))
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "contig_vs_n50.scatter.pdf"),
  plot = ptScat_contig_vs_n50, width = 8, height = 8
)

## #contigs vs L50
ptScat_contig_vs_l50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = L50, y = n_contigs)
) +
  geom_point() +
  labs(title = "#contigs vs L50") +
  pt_theme

ggsave(
  filename = file.path(outDir, "contig_vs_l50.scatter.pdf"),
  plot = ptScat_contig_vs_l50, width = 8, height = 8
)


## L50 vs N50
ptScat_l50_vs_n50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = N50, y = L50)
) +
  geom_point() +
  labs(title = "L50 vs N50") +
  scale_x_continuous(
    labels = label_comma(scale_cut = c(bp = 0, kb = 1000, mb = 1000000))
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "l50_vs_n50.scatter.pdf"),
  plot = ptScat_l50_vs_n50, width = 8, height = 8
)



#####################################################################
# busco


ptHist_buscoSC <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = complete_single_copy)
) +
  geom_histogram(bins = 100, color = "black", fill = "black") + 
  labs(
    title = "Histogram of Single Copy BUSCO (excluding 2 outliers)",
    x = "Single Copy BUSCO", y = "Count"
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "BUSCO_SC.hist.pdf"),
  plot = ptHist_buscoSC, width = 8, height = 8
)


## single copy BUSCO vs N50
ptScat_buscoSC_vs_n50 <- ggplot2::ggplot(
  data = dplyr::filter(.data = genomeMetadata, n_contigs < 500),
  mapping = aes(x = N50, y = complete_single_copy)
) +
  geom_point() +
  labs(title = "single copy BUSCO vs N50") +
  scale_x_continuous(
    labels = label_comma(scale_cut = c(bp = 0, kb = 1000, mb = 1000000))
  ) +
  pt_theme

ggsave(
  filename = file.path(outDir, "BUSCO_SC_vs_n50.scatter.pdf"),
  plot = ptScat_buscoSC_vs_n50, width = 8, height = 8
)


