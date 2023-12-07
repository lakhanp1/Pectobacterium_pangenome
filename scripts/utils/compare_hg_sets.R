
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
  
  fwSynLen <- 0
  revSynLen <- 0
  if (!is.null(forwardSyn)) {
    fwSynLen <- length(forwardSyn$lcs)
  }
  
  if (!is.null(reverseSyn)) {
    revSynLen <- length(reverseSyn$lcs)
  }
  
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
  lcs <- NULL
  lcsFound <- FALSE
  
  thisLcs <- c()
  lcsAln <- list(s1 = vector(mode = class(seqA)), s2 = vector(mode = class(seqB)))
  lcsPos <- list(s1 = numeric(0L), s2 = numeric(0L))
  cigar <- c()
  
  i <- m + 1
  j <- n + 1
  gapLen <- 0
  lcsScore <- 0
  highestScore <- 0
  
  while (i > 1 && j > 1) {
    
    # if(seqA[i-1] == seqB[j-1]){
    if(tb[i, j] == 0){
      thisLcs <- c(seqA[i-1], thisLcs)
      lcsAln$s1 <- c(seqA[i-1], lcsAln$s1)
      lcsAln$s2 <- c(seqB[j-1], lcsAln$s2)
      lcsPos$s1 <- c(i-1, lcsPos$s1)
      lcsPos$s2 <- c(j-1, lcsPos$s2)
      
      gapLen <- 0
      if(highestScore < dp[i, j]){
        highestScore <- dp[i, j]
      }
      
      lcsScore <- highestScore - dp[i-1, j-1]
      
      i <- i - 1
      j <- j - 1
      
      # } else if(dp[i, j] == (dp[i-1, j] + gapPenalty)){
    } else if(tb[i, j] == -1){
      lcsAln$s1 <- c(seqA[i-1], lcsAln$s1)
      lcsAln$s2 <- c("-", lcsAln$s2)
      gapLen <- gapLen + 1
      
      i <- i - 1
      
    } else if(tb[i, j] == 1){
      lcsAln$s1 <- c("-", lcsAln$s1)
      lcsAln$s2 <- c(seqB[j-1], lcsAln$s2)
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
          lcsFound <- TRUE
        }
        
        # optionally, save a valid lcs in list
        lcsList <- append(
          lcsList,
          list(list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore))
        )
      }
      thisLcs <- c()
      lcsAln <-  list(s1 = vector(mode = class(seqA)), s2 = vector(mode = class(seqB)))
      lcsPos <- list(s1 = numeric(0L), s2 = numeric(0L))
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
        lcsAln <-  list(s1 = vector(mode = class(seqA)), s2 = vector(mode = class(seqB)))
        lcsPos <- list(s1 = numeric(0L), s2 = numeric(0L))
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
      lcsFound <- TRUE
    }
    
    # optionally, save a valid lcs in list
    lcsList <- append(
      lcsList,
      list(list(lcs = thisLcs, aln = lcsAln, pos = lcsPos, score = lcsScore))
    )
    
  }
  
  # include the CIGAR string in LCS results
  # M: alignment match (can be a sequence match or mismatch)
  # I: insertion to the reference
  # D: deletion from the reference
  # N: skipped region from the reference
  # S: soft clipping (clipped sequences present in SEQ)
  # H: hard clipping (clipped sequencesT present in SEQ)
  # P: padding (silent deletion from padded reference)
  # =: sequence match
  # X: sequence mismatch
  
  
  # seqA
  # seqB
  # lcsList
  # purrr::map_chr(lcsList$aln, paste, collapse = "|")
  
  if(lcsFound){
    return(lcs)
  } else {
    return(NULL)
  }
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
