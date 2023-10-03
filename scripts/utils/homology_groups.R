#' extract homology groups for set of genomes in pangenome
#'
#' @param file pangenome homology group file
#' @param genomes genomes to select for
#' @param groups group category to return. One of `c("core", "accessory",
#'  "unique", "all")`. default: "all"
#' @param pav Whether to convert homology group count >1 to 1 i.e. output
#' binary PAV matrix. Default: FALSE
#'
#' @return A homology group tibble
#' @export
#'
#' @examples NA
homology_groups_extract <- function(file, genomes, groups = "all", pav = FALSE) {
  stopifnot(
    groups %in% c("core", "accessory", "unique", "all")
  )
  
  ## homology groups PAV matrix
  hgs <- suppressMessages(readr::read_csv(file)) %>%
    dplyr::rename_with(.fn = ~ stringr::str_replace_all(.x, "( |-)", "_")) %>%
    dplyr::rename_with(.fn = ~ stringr::str_replace_all(.x, "Genome_", "")) %>%
    dplyr::rename_with(.fn = ~ tolower(.x)) %>%
    dplyr::mutate(homology_group_id = as.character(homology_group_id)) %>%
    dplyr::select(hg = homology_group_id, tidyselect::all_of(genomes))
  
  ## binary matrix for homology_group PAV
  hgBinaryMat <- dplyr::mutate(
    hgs,
    dplyr::across(
      .cols = tidyselect::all_of(genomes),
      .fns = ~ dplyr::if_else(.x > 1, true = 1, false = .x)
    )
  ) %>%
    tibble::column_to_rownames(var = "hg") %>%
    as.matrix()
  
  hgSum <- matrixStats::rowSums2(
    x = hgBinaryMat, useNames = T
  ) %>%
    tibble::enframe(name = "hg", value = "nGenomes") %>%
    dplyr::filter(nGenomes != 0) %>%
    dplyr::mutate(
      class = dplyr::case_when(
        nGenomes == 1 ~ "unique",
        nGenomes == !!nrow(sampleInfo) ~ "core",
        nGenomes < !!nrow(sampleInfo) & nGenomes > 1 ~ "accessory"
      )
    )
  
  if (groups != "all") {
    hgSum %<>% dplyr::filter(class == .env$groups)
  }
  
  if (pav) {
    hgs <- tibble::as_tibble(hgBinaryMat, rownames = "hg")
  }
  
  hgSet <- dplyr::left_join(hgSum, hgs, by = "hg") %>%
    dplyr::select(-nGenomes)
  
  return(hgSet)
}

################################################################################

#' prepare homology group PAV matrix from pan.db
#'
#' @param pandb A OrgDb object that store pangenome data
#' @param type One of c("pav", "cnv")
#' @param groups Homology groups to return. Default: all groups
#'
#' @return a matrix of homology groups PAV/CNV per genome
#' @export
#'
#' @examples
homology_groups_mat <- function(pandb, type, groups = NULL) {
  type <- match.arg(arg = tolower(type), choices = c("pav", "cnv"))
  
  stopifnot(
    isa(pandb, c("OrgDb"))
  )
  
  hgs <- AnnotationDbi::select(
    x = pandb, keys = keys(pandb), columns = c("genomeId", "mRNA_key")
  )
  
  if (type == "cnv") {
    hgWide <- tidyr::pivot_wider(
      hgs,
      names_from = "GID", values_from = "mRNA_key",
      values_fn = ~ length(.x), values_fill = 0
    )
  } else {
    hgWide <- tidyr::pivot_wider(
      hgs,
      names_from = "GID", values_from = "mRNA_key",
      values_fn = ~ length(.x[1]), values_fill = 0
    )
  }
  
  if (!is.null(groups)) {
    hgWide %<>% dplyr::select(genomeId, tidyselect::all_of(groups))
  }
  
  hgMat <- tibble::column_to_rownames(hgWide, var = "genomeId") %>%
    as.matrix()
  
  return(hgMat)
}

################################################################################
#' Plot homology group heatmap with provided clustering
#'
#' @param mat homology group matrix
#' @param phy A phylo object
#' @param metadata metadata for genomes. Must have SpeciesName and Genome column
#' @param speciesInfo Optionally, a metadata for genomes. Must have SpeciesName and
#' Genome column. In addition to the genome similarity heatmap, a species key is
#' shown as heatmap in the pangenome
#' @param hgAn a data.frame for homology group annotation. It should have three
#' columns: `c(hg_id, hg_group, hg_label)`. Default: `NULL`
#' @param markGenomes A list with two elements named `list(compare = c(),
#' against = c())`. Default: `NULL`
#'
#' @return ComplexHeatmap with layout: species key heatmap | ANI heatmap
#' @export
#'
#' @examples
homology_group_heatmap <- function(mat, phy, speciesInfo = NULL,
                                   hgAn = NULL, markGenomes = NULL, ...) {
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    # https://github.com/jokergoo/ComplexHeatmap/issues/949
    all(rownames(mat) == phy$tip.label),
    is.null(speciesInfo) | (all(rownames(mat) %in% speciesInfo$genomeId)),
    is.null(hgAn) | all(colnames(mat) == hgAn$hg_id),
    any(class(phy) == "phylo"),
    is.null(markGenomes) |
      (is.list(markGenomes) & all(unlist(markGenomes) %in% rownames(mat)))
  )
  
  ## homology groups heatmap arguments
  ht_args <- list(
    matrix = mat,
    col = structure(
      viridisLite::viridis(n = max(3, min(max(mat), 6)) + 1, option = "B"),
      names = seq(0, max(3, min(max(mat), 7)))
    ),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    show_row_names = FALSE,
    cluster_columns = TRUE, cluster_column_slices = FALSE,
    show_column_names = TRUE, column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 10),
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    ...
  )
  
  if (!is.null(hgAn)) {
    ht_args$column_split <- hgAn$hg_group
    ht_args$column_order <- hgAn$hg_id
    ht_args$column_labels <- hgAn$label
    ht_args$cluster_columns <- FALSE
  }
  
  ## homology group heatmap
  ht_hg <- do.call(
    what = ComplexHeatmap::Heatmap,
    args = ht_args
  )
  
  # draw(ht_hg)
  htList <- ht_hg
  
  # optional species key heatmap
  if (!is.null(speciesInfo)) {
    stopifnot(
      all(rownames(mat) %in% speciesInfo$genomeId),
      has_name(speciesInfo, "genomeId"),
      has_name(speciesInfo, "SpeciesName")
    )
    
    ht_species <- species_key_heatmap(
      genomes = rownames(mat), speciesInfo = speciesInfo,
      markGenomes = markGenomes
    )
    
    htList <- ht_species + ht_hg
  }
  
  return(htList)
}

################################################################################

#' Get homology groups for a genomic region
#'
#' @param pandb org.db pangenome object
#' @param genome genome identifier in pangenome
#' @param chr chromosome name
#' @param start start. default: `1` i.e. start of chromosome
#' @param end end position. default: `Inf` end of chromosome
#'
#' @return A vector of homology group identifiers
#' @export
#'
#' @examples
region_homology_groups <- function(pandb, genome, chr, start = 1, end = Inf) {
  stopifnot(
    !is.na(chr)
  )
  
  start <- ifelse(is.na(start), 1, start)
  end <- ifelse(is.na(end), Inf, end)
  
  df <- suppressMessages(AnnotationDbi::select(
    x = pandb, keys = genome, keytype = "genomeId",
    columns = c("GID", "chr_name", "start", "end")
  )) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(start, end), .fns = as.integer
      )
    ) %>%
    dplyr::filter(chr_name == !!chr, start >= !!start, end <= !!end) %>%
    dplyr::arrange(start)
  
  return(df$GID)
}
################################################################################

#' Locate homology group set position on a contig or chromosome
#'
#' @param hgs homology groups
#' @param genome Genome against which to map
#' @param pandb org.db pangenome object
#' @param chr chromosome or contig name to restrict
#'
#' @return A number `[0, 1]` giving relative position of homology groups on contig.
#' `-1` if the homology groups are scattered on a contig.
#' @export
#'
#' @examples
get_hg_sets_location <- function(hgs, genome, chr, pandb) {
  hgInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = hgs,
      columns = c("genomeId", "chr", "chr_id", "chr_name", "start", "end", "strand")
    )
  ) %>%
    dplyr::filter(genomeId == !!genome, chr_name == !!chr) %>%
    dplyr::mutate(dplyr::across(c(start, end), as.integer)) %>%
    dplyr::arrange(start)
  
  stopifnot(length(unique(hgInfo$chr_id)) == 1)
  
  chrInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = genome, keytype = "genomeId",
      columns = c("GID", "genomeId", "chr", "chr_id", "chr_name", "start", "end", "strand")
    )
  ) %>%
    dplyr::filter(chr_id == !!hgInfo$chr_id[1]) %>%
    dplyr::mutate(dplyr::across(c(start, end), as.integer)) %>%
    dplyr::arrange(start)
  
  # use RLE structure to check for tandem overlap
  hgRle <- rle(chrInfo$GID %in% hgInfo$GID)
  
  if (!any(hgRle$lengths[hgRle$values == TRUE] >= nrow(hgInfo))) {
    # for a tandem match of homology groups, there should be one (RLE length for TRUE) >= #HGs
    # if these is no such RLE length, something is wrong
    return(-1)
  } else {
    if (identical(hgRle$values, c(TRUE))) {
      # hgs are on independent contig
      return(0)
    } else if (identical(hgRle$values[1:2], c(TRUE, FALSE)) &
               hgRle$lengths[1] >= nrow(hgInfo)) {
      # hgs are at the beginning
      return(0)
    } else if (identical(tail(hgRle$values, 2), c(FALSE, TRUE)) &
               tail(hgRle$lengths, 1) >= nrow(hgInfo)) {
      # hgs are at the beginning
      return(1)
    } else {
      # hgs present within a contig
      
      hgPosition <- which(
        hgRle$values == TRUE & hgRle$lengths >= nrow(hgInfo)
      )
      
      return((sum(head(hgRle$lengthsm, hgPosition - 1)) + 1) / nrow(chrInfo))
    }
  }
}


################################################################################


#' From a set of homology groups, detect homology groups that are arranged in
#' tandem in a genome
#'
#' @param hgs A vector of homology group ids
#' @param pandb org.db pangenome object
#' @param genomeId Find tandem match against a specific genome
#'
#' @return A list of homology groups where each element is a tandem group set
#' @export
#'
#' @examples
tandem_hg_match <- function(hgs, pandb, genomeId = NULL) {
  hgInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = hgs,
      columns = c("genomeId", "chr", "chr_id", "chr_name", "start", "end", "strand")
    )
  ) %>%
    dplyr::mutate(dplyr::across(c(start, end), as.integer)) %>%
    dplyr::arrange(start)
}


################################################################################

#' Check if two homology group sets overlap with synteny conserved
#'
#' Implement https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1626491/ algorithm
#'
#' @param ref reference homology group ids order as per genomic coordinates
#' @param qur query homology group ids order as per genomic coordinates
#' @param ... other arguments to `longest_local_subsequence()`
#' 
#' @return
#' @export
#'
#' @examples NA
syntenic_hg_overlap <- function(ref, qur, ...) {
  
  forwardSyn <- longest_local_subsequence(seq1 = ref, seq2 = qur, ...)
  reverseSyn <- longest_local_subsequence(seq1 = ref, seq2 = rev(qur), ...)
  
  # fwSynLen <- purrr::map_dbl(forwardSyn, ~length(.x$lcs)) %>% sum()
  # revSynLen <- purrr::map_dbl(reverseSyn, ~length(.x$lcs)) %>% sum()
  fwSynLen <- length(forwardSyn$lcs)
  revSynLen <- length(reverseSyn$lcs)
  
  if(fwSynLen > revSynLen){
    return(forwardSyn)
  } else{
    return(reverseSyn)
  }
  
}

################################################################################

#' Find a longest common subsequence between the two arrays
#'
#' @param seq1 array 1
#' @param seq2 array 2
#' @param matchScore match score. Default: 5
#' @param gapPenalty gap open penalty. Default: -2
#' @param maxGapLen longest possible gap. Default: 2
#' @param minChainLen minimum tandem syntenic matches. Default: 5
#'
#' @return a list with two elements: longest local subsequence between two arrays
#' and position of this subsequence in seq1
#' 
#' @export
#'
#' @examples 
#' seq1 <- c(1, 2, 3, 4, 5, 6)
#' seq2 <- c(9, 3, 4, 5, 8)
#' lcs <- longest_local_subsequence(seq1, seq2)
#' cat("Longest local subsequence:", lcs$lcs)
#' 
longest_local_subsequence <- function(
    seq1, seq2, matchScore = 5, gapPenalty = -2, maxGapLen = 2, minChainLen = 5
) {
  # matchScore            -Z: a constant match score. can also use E-value
  # gapPenalty            -o (GAP_OPEN_PENALTY) in DAGchainer 
  # maxGapLen             -D (MAX_DIST_BETWEEN_MATCHES) in DAGchainer
  # minChainLen           -A in DAGchainer
  
  # gapExtPenalty <- -2   # -e (INDEL_SCORE) in DAGchainer 
  # minAlnScore <- matchScore * minChainLen   # -x (MIN_ALIGNMENT_SCORE) in DAGchainer
  
  # seq1 <- c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "r", "s", "t")
  # seq2 <- c("p", "q", "a", "b", "c", "d", "t", "f", "g", "h", "i", "l", "m", "n", "x", "y", "z", "i", "j", "k")
  
  # use the longest sequence as seqA
  seqA <- seq1
  seqB <- seq2
  # if(length(seq1) < length(seq2)){
  #   seqA <- seq2
  #   seqB <- seq1
  # }
  
  m <- length(seqA)
  n <- length(seqB)
  
  # empty scoring matrix and traceback matrix
  dp <- matrix(0, nrow = m + 1, ncol = n + 1, dimnames = list(c(".", seqA), c(".", seqB)))
  tb <- matrix(0, nrow = m + 1, ncol = n + 1, dimnames = list(c(".", seqA), c(".", seqB)))
  maxDpScore <- 0
  maxScoreRow <- 0
  maxScoreCol <- 0
  
  # populate the matrix using dynamic programming
  for (i in 1:(m+1)) {
    for (j in 1:(n+1)) {
      
      if(i > 1 && j > 1){
        if(seqA[i-1] == seqB[j-1]){
          
          dp[i, j] <- dp[i-1, j-1] + matchScore
          tb[i, j] <- 0
          
        } else{
          
          dp[i, j] <- max(
            dp[i-1, j] + gapPenalty,
            dp[i, j-1] + gapPenalty
          )
          
          if(dp[i, j] == dp[i-1, j] + gapPenalty){
            tb[i, j] <- -1
          } else{
            tb[i, j] <- 1
          }
          
        }
      }
      
      if(dp[i, j] > maxDpScore){
        maxDpScore <- dp[i, j]
        maxScoreRow <- i
        maxScoreCol <- j
      }
    }
  }
  
  # dp
  
  # traceback to reconstruct the longest local subsequence
  lcsList <- list()
  lcs <- list(lcs = c()) 
  
  thisLcs <- c()
  lcsAln <- list(vector(mode = class(seqA)), vector(mode = class(seqB)))
  lcsPos <- list(numeric(0L), numeric(0L))
  
  i <- m + 1
  j <- n + 1
  gapLen <- 0
  lcsScore <- 0
  highestScore <- 0
  
  while (i > 1 && j > 1) {
    
    # if(seqA[i-1] == seqB[j-1]){
    if(tb[i, j] == 0){
      thisLcs <- c(seqA[i-1], thisLcs)
      lcsAln[[1]] <- c(seqA[i-1], lcsAln[[1]])
      lcsAln[[2]] <- c(seqB[j-1], lcsAln[[2]])
      lcsPos[[1]] <- c(i-1, lcsPos[[1]])
      lcsPos[[2]] <- c(j-1, lcsPos[[2]])
      
      gapLen <- 0
      if(highestScore < dp[i, j]){
        highestScore <- dp[i, j]
      }
      
      lcsScore <- highestScore - dp[i-1, j-1]
      
      i <- i - 1
      j <- j - 1
      
      # } else if(dp[i, j] == (dp[i-1, j] + gapPenalty)){
    } else if(tb[i, j] == -1){
      lcsAln[[1]] <- c(seqA[i-1], lcsAln[[1]])
      lcsAln[[2]] <- c("-", lcsAln[[2]])
      gapLen <- gapLen + 1
      
      i <- i - 1
      
    } else if(tb[i, j] == 1){
      lcsAln[[1]] <- c("-", lcsAln[[1]])
      lcsAln[[2]] <- c(seqB[j-1], lcsAln[[2]])
      gapLen <- gapLen + 1
      
      j <- j - 1
      
    }
    
    # add thisLcs to lcsList if the gap > maxGapLen
    if(gapLen > maxGapLen){
      # remove the gaps at the begining
      lcsAln <- purrr::map(.x = lcsAln, .f = ~ .x[-c(1:gapLen)])
      
      if(length(thisLcs) >= minChainLen){
        # update the longest LCS
        if(length(thisLcs) > length(lcs$lcs)){
          lcs <- list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore)
        }
        
        # optionally, save a valid lcs in list
        lcsList <- append(
          lcsList,
          list(list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore))
        )
      }
      thisLcs <- c()
      lcsAln <- list(vector(mode = class(seqA)), vector(mode = class(seqB)))
      lcsPos <- list(numeric(0L), numeric(0L))
      lcsScore <- 0
      highestScore <- 0
      
    }
    
    # cat("i=", i, "; j=", j, "\n")
    # if traceback never went to the highest scoring address, 
    # set i and j to explicitly traceback from the highest scoring address
    if(j == 1 || i == 1){
      if(i > maxScoreRow || j > maxScoreCol){
        i <- maxScoreRow
        j <- maxScoreCol
        
        thisLcs <- c()
        lcsAln <- list(vector(mode = class(seqA)), vector(mode = class(seqB)))
        lcsPos <- list(numeric(0L), numeric(0L))
        lcsScore <- 0
        highestScore <- 0
      }
    }

    
  }
  
  # add final LCS to lcsList
  if(length(thisLcs) >= minChainLen){
    # update the longest LCS
    if(length(thisLcs) > length(lcs$lcs)){
      lcs <- list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore)
    }
    
    # optionally, save a valid lcs in list
    lcsList <- append(
      lcsList,
      list(list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore))
    )
    
  }
  
  
  # seqA
  # seqB
  # lcsList
  # purrr::map_chr(lcsList$aln, paste, collapse = "|")
  
  return(lcs)
}

################################################################################

#' Print the `longest_local_subsequence()` results
#'
#' @param lcs `longest_local_subsequence()` results
#'
#' @return NULL
#' @export
#'
#' @examples NA
print_lcs <- function(lcs){

    cat(
      "LCS: ", paste(lcs$lcs, collapse = " "),
      "\nLength: ", length(lcs$lcs),
      "\nScore: ", lcs$score, 
      "\n", sep = ""
    )
    
    cat("Alignment:\n")
    purrr::map_chr(lcs$aln, paste, collapse = " ") %>% cat(sep = "\n")
    cat("\n")

  

}

################################################################################
