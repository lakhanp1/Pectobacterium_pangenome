suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))

# process and combine genomad and checkv output files for all samples

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
################################################################################

confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]

################################################################################

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId 
)

# sampleInfo %<>% dplyr::filter(sampleId %in% c("GCF_004296765.1_ASM429676v1", "JGAR-100719-1"))
# sampleId <- "GCF_009931535.1_ASM993153v1"
prophageDf <- NULL
plasmidDf <- NULL

for (sampleId in sampleInfo$sampleId) {
  genomadPath <- paste(confs$data$prophages$dir, "/", sampleId, "/", sampleId, sep = "")
  
  cat(genomadPath, "...\n")
  
  ## prophage data
  virSummary <- NULL
  
  virSummary <- suppressMessages(
    readr::read_tsv(
      file = paste(genomadPath, "_summary/", sampleId, "_virus_summary.tsv", sep = "")
    )
  ) %>%
    dplyr::select(-n_genes) %>%
    dplyr::rename_with(
      .fn = ~ paste("genomad.", .x, sep = ""),
      .cols = c(virus_score, fdr, n_hallmarks, marker_enrichment)
    ) %>%
    dplyr::mutate(
      sampleId = !!sampleId,
      genomeId = !!sampleInfoList[[sampleId]]$genomeId,
      SpeciesName = !!sampleInfoList[[sampleId]]$SpeciesName,
      chr = stringr::str_replace(seq_name, "(.*)\\|.*", "\\1"),
      .before = coordinates
    ) %>%
    tidyr::separate(
      col = coordinates, into = c("start", "end"), sep = "-", convert = TRUE
    )
  
  if (nrow(virSummary) > 0) {
    
    # import checkv contamination and summary
    vContamination <- suppressMessages(
      readr::read_tsv(
        file = paste(genomadPath, "_checkv/contamination.tsv", sep = "")
      )
    ) %>% 
      dplyr::select(
        contig_id, contamination_types = region_types,
        contamination_regions = region_coords_bp
      ) %>% 
      dplyr::filter(!is.na(contamination_types)) %>% 
      dplyr::mutate(host_contamination = 1)
    
    if (any(!vContamination$contamination_types %in% c("host,viral", "viral,host", "host,viral,host"))) {
      stop(
        "Unusual virus contamination type in ", sampleId, ": ",
        paste(vContamination$contamination_types, collapse = "; ")
      )
    }
    
    hostFilteredRegion <- dplyr::mutate(
      vContamination,
      dplyr::across(
        .cols = c(contamination_types, contamination_regions),
        .fns = ~stringr::str_split(.x, pattern = ",")
      )
    ) %>% 
      tidyr::unnest(cols = c(contamination_types, contamination_regions)) %>% 
      dplyr::filter(contamination_types == "viral") %>% 
      tidyr::separate(
        col = contamination_regions, into = c("vStart", "vEnd"),
        sep = "-", convert = TRUE
      ) %>% 
      dplyr::select(contig_id, vStart, vEnd)
    
    vContamination %<>% dplyr::left_join(
      y = hostFilteredRegion, by = "contig_id"
    )
    
    checkv <- suppressMessages(
      readr::read_tsv(
        file = paste(genomadPath, "_checkv/quality_summary.tsv", sep = "")
      )
    ) %>%
      dplyr::left_join(vContamination, by = "contig_id") %>% 
      tidyr::replace_na(replace = list(host_contamination = 0)) %>% 
      dplyr::select(-proviral_length)
    
    # update the prophage coordinates if there is a contamination
    virSummary <- dplyr::left_join(
      virSummary, checkv, by = c("seq_name" = "contig_id")
    ) %>% 
      dplyr::mutate(
        across(.cols = c(start, end), .fns = ~ as.double(.x)),
        start = if_else(
          condition = host_contamination == 1,
          true = start + vStart - 1, false = start
        ),
        end = if_else(
          condition = host_contamination == 1,
          true = start + vEnd - vStart, false = end
        ),
        prophage_length = end - start + 1
      ) %>% 
      dplyr::mutate(
        prophage_id = paste(genomeId, ".vir_", 1:n(), sep = ""),
        prophage_length = dplyr::if_else(
          condition = is.na(prophage_length), length, prophage_length
        )
      ) %>% 
      dplyr::select(-length)
    
    prophageDf <- dplyr::bind_rows(prophageDf, virSummary)
    
  }
  
  ## plasmid data
  plasmidSummary <- NULL
  
  plasmidSummary <- suppressMessages(
    readr::read_tsv(
      file = paste(genomadPath, "_summary/", sampleId, "_plasmid_summary.tsv", sep = "")
    )
  ) %>% 
    dplyr::rename(chr = seq_name, plasmid_length = length) %>% 
    dplyr::mutate(
      sampleId = !!sampleId,
      genomeId = !!sampleInfoList[[sampleId]]$genomeId,
      SpeciesName = !!sampleInfoList[[sampleId]]$SpeciesName,
      .before = chr
    )
  
  if (nrow(plasmidSummary) > 0) {
    plasmidSummary %<>%  dplyr::mutate(
      plasmid_id = paste(genomeId, ".pmd_", 1:n(), sep = "")
    )
    
    plasmidDf <- dplyr::bind_rows(plasmidDf, plasmidSummary)
  }
  
}

## save prophage data
prophageDf <- dplyr::select(
  prophageDf, contig_id = seq_name, prophage_id, sampleId, SpeciesName, chr, start, end,
  prophage_length, topology, taxonomy, starts_with("genomad."), everything(),
  -vStart, -vEnd
)

readr::write_tsv(prophageDf, file = confs$data$prophages$files$data)

dplyr::select(prophageDf, sampleId, chr, start, end, prophage_id) %>%
  dplyr::mutate(
    fna = paste(confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""),
    region = dplyr::if_else(
      is.na(start), true = chr, false = paste(chr, ":", start, "-", end, sep = ""),
    ),
    out = paste(confs$data$prophages$dir, "/phage_seqs/", prophage_id, ".fna", sep = ""),
    faidx = paste("samtools faidx", fna, region, ">", out)
  ) %>%
  dplyr::select(faidx) %>%
  readr::write_tsv(
    file = "./scripts/preprocessing/prophage_genomes_extract.sh",
    col_names = FALSE
  )

plasmidDf <- dplyr::select(
  plasmidDf, plasmid_id, everything()
)

## save plasmid data
readr::write_tsv(plasmidDf, file = confs$data$plasmids$files$data)

dplyr::select(plasmidDf, sampleId, chr, plasmid_id) %>%
  dplyr::mutate(
    fna = paste(confs$data$prokka$dir, "/", sampleId, "/", sampleId, ".fna", sep = ""),
    region = chr,
    out = paste(confs$data$plasmids$dir, "/plasmid_seqs/", plasmid_id, ".fna", sep = ""),
    faidx = paste("samtools faidx", fna, region, ">", out)
  ) %>%
  dplyr::select(faidx) %>%
  readr::write_tsv(
    file = "./scripts/preprocessing/plasmid_genomes_extract.sh",
    col_names = FALSE
  )


