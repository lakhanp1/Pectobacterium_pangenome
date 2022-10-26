suppressMessages(library(tidyverse))
suppressMessages(library(magrittr))
suppressMessages(library(rtracklayer))
suppressMessages(library(Biostrings))
suppressMessages(library(GenomicRanges))
suppressMessages(library(BSgenome))



# Author: Lakhansing Pardeshi
# Contact: lakhansing.pardeshi@wur.nl
## Find duplicate genomes assemblies and detect positions where there are mismatches

rm(list = ls())

#####################################################################
file_chrSize <- here::here("analysis", "01_multiqc", "assembly_chr_size.txt")
file_summary <- here::here("data/reference_data", "sample_metadata.tsv")
path_2bit <- here::here("temp", "local_copy")

#####################################################################

filteredAssemblies <- suppressMessages(readr::read_tsv(file = file_summary)) %>% 
  dplyr::select(
    sampleId, ExclFromRefSeq, Anomalous, replaced
  ) %>%
  dplyr::filter(
    !is.na(ExclFromRefSeq) | !is.na(Anomalous) | !is.na(replaced)
  ) 


chrSizes <- suppressMessages(readr::read_tsv(
  file = file_chrSize, col_names = c("sampleId", "chr", "length"))
)


chrL <- dplyr::group_by(chrSizes, sampleId) %>% 
  dplyr::arrange(desc(length), .by_group = TRUE) %>% 
  tidyr::nest() %>% 
  # head() %>% 
  purrr::transpose() %>% 
  purrr::map(
    .f = function(a){
      list(
        sampleId = a$sampleId,
        length = a$data$length
      )
    }
  )

chrL %<>% purrr::discard(.p = function(x) x$sampleId %in% filteredAssemblies$sampleId)


## function to compare two fasta files where chromosome lengths are exactly same
verify_duplicate_assemblies <- function(s1, s2){
  sq1 = rtracklayer::import.2bit(con = s1)
  sq2 = rtracklayer::import.2bit(con = s2)
  
  # dna1 <- DNAStringSet(c("AAA", "TC", "CCAT", "TAC", "AAA", "CAAC", "CG"))
  # dna2 <- DNAStringSet(c("AAA", "TC", "CATG", "TACA", "AAA", "CAC", "GC"))
  # match(dna1, dna2)
  # sort(dna1) %in% sort(dna2)
  
  if(all(sort(sq1) %in% sort(sq2))){
    
    return(tibble::tibble(identical = TRUE))
    
  } else{
    
    ## split big chromosomes into 10kb GRs
    sq1Bins <- GenomicRanges::tileGenome(
      seqlengths = seqinfo(sq1), tilewidth = 10000, cut.last.tile.in.chrom = T
    )
    
    sq2Bins <- GenomicRanges::tileGenome(
      seqlengths = seqinfo(sq2), tilewidth = 10000, cut.last.tile.in.chrom = T
    )
    
    ## get respective sequences for chromosomes
    sq1BinSeq <- BSgenome::getSeq(x = sq1, names = sq1Bins)
    sq2BinSeq <- BSgenome::getSeq(x = sq2, names = sq2Bins)
    
    mismatchedBins <- which(sq1BinSeq != sq2BinSeq)
    
    ## pairwise alignment and mismatch summary table
    pa <- Biostrings::pairwiseAlignment(
      pattern = sq1BinSeq[mismatchedBins], subject = sq2BinSeq[mismatchedBins]
    )
    
    ## select the bins which did not match to detect mismatched positions
    mismatches <- Biostrings::mismatchTable(x = pa)
    
    ## change the start and end coordinates relative to the complete chromosome
    mismatches %<>%
      dplyr::mutate(
        identical = FALSE,
        n_mismatch = n(),
        dplyr::across(
          .cols = c(ends_with("Start"), ends_with("End")),
          .fns = ~ .x + GenomicRanges::start(sq1Bins[mismatchedBins])
        ),
        PatternChr = as.character(GenomicRanges::seqnames(sq1Bins[mismatchedBins])),
        SubjectChr = as.character(GenomicRanges::seqnames(sq2Bins[mismatchedBins])),
      ) %>% 
      dplyr::select(
        identical, n_mismatch, starts_with("Pattern"), starts_with("Subject"),
        -PatternId
      )
    
    return(mismatches)
  }
  
}


duplicateGenomes <- tibble::tibble(genome1 = character(0), genome2 = character(0))


for (i in 1:length(chrL)) {
  for (j in 1:i) {
    if(i == j) next
    
    if(length(chrL[[i]]$length) == length(chrL[[j]]$length)){
      if(all(chrL[[i]]$length == chrL[[j]]$length)){
        
        sample1 = chrL[[i]]$sampleId
        sample2 = chrL[[j]]$sampleId
        
        duplicateGenomes %<>% dplyr::bind_rows(
          tibble::tibble(
            genome1 = sample1,
            genome2 = sample2,
            duplicate = verify_duplicate_assemblies(
              s1 = paste(path_2bit, "/", sample1, ".2bit", sep = ""),
              s2 = paste(path_2bit, "/", sample2, ".2bit", sep = "")
            )
          )
        )
      }
    }
  }
}

duplicateGenomes <- tidyr::unnest(duplicateGenomes, cols = c(duplicate))

readr::write_tsv(
  x =duplicateGenomes,
  file = here::here("analysis", "01_multiqc", "duplicate_genomes.tab")
)



