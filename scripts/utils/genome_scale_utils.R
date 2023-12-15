
#' get orderd genome identifiers based on genome similarity
#'
#' @param genome genome identifier in the pangeonme
#' @param orgdb pangenome org.db object
#' @param n how many closest genomes to return. Default: `NULL` i.e. return all
#' genomes ranked based on closeness
#'
#' @return A vector
#' @export
#'
#' @examples
genome_similarity <- function(genome, orgdb, n = NULL){
  
}


################################################################################

#' Invert the coordinates for a region
#' 
#' For the circular genomes (bacteria, virus, plasmids, etc.), the genome assembly
#' orientation can differ for the same species genomes. Comparing or visualizing
#' regions from such genomes can give false sense of inversion. Therefore, this
#' function convert the coordinates of such genomic regions provided as a dataframe.
#'
#' @param df A genomic region dataframe with `start`, `end` and `strand` columns.
#'
#' @return Inverted coordinates dataframe
#' @export
#'
#' @examples
#' x <- tibble::tibble(
#'   start = c(2, 9, 16, 21, 33),
#'   end = c(6, 11, 19, 29, 40),
#'   strand = c(1, -1, 1, 1, -1),
#'   name = LETTERS[1:5]
#' )
#' 
#' invert_coordinates(x)
#' 
invert_coordinates <- function(df){
  stopifnot(
    hasName(df, c("start", "end", "strand")),
    is.numeric(df$start),
    is.numeric(df$end),
    is.numeric(df$strand)
  )
  
  # x <- tibble::tibble(
  #     start = c(2, 9, 16, 21, 33),
  #     end = c(6, 11, 19, 29, 40),
  #     strand = c(1, -1, 1, 1, -1),
  #     name = LETTERS[1:5]
  #   )

  df <- dplyr::arrange(df, start) %>% 
    dplyr::mutate(
      "_strand" = strand,
      "_start" = start,
      "_end" = end,
      strand = strand * -1,
      width = end - start
    )
  
  df$start <- diffinv(
    x = (diff(df$end) * -1),
    xi = tail(df$end, n = 1) - df$width[1]
  )
  
  df$end <- df$start + df$width
  df$width <- NULL
  
  df <- dplyr::arrange(df, start)
  return(df)
}

################################################################################

