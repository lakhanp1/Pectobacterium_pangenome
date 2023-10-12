
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

#' import prophage relationships as igraph graph
#'
#' @param file prophage relationship file
#' @param nodeAn node annotations
#' @param key column to use for combining node annotation
#'
#' @return an `igraph` graph object
#' @export
#'
#' @examples NA
prophage_igraph <- function(file, nodeAn, key){
  
  phageRelations <- suppressMessages(readr::read_tsv(file))
  
  stopifnot(
    key %in% colnames(phageRelations),
    'parent' %in% colnames(phageRelations),
    key %in% colnames(nodeAn)
  )
  
  nodes <- dplyr::select(phageRelations, !!key, nodeType, genomeId,
                         nHgs, jaccardIndex) %>% 
    dplyr::left_join(y = nodeAn, by = key) %>% 
    dplyr::rename(id = !!key)
  
  edges <- dplyr::filter(phageRelations, !is.na(parent)) %>% 
    dplyr::select(
      from = !!key, to = parent, weight = jaccardIndex,
      nSharedHgs, nSyntenicSharedHgs, starts_with("perShared"), relation
    )
  
  
  if(length(setdiff(union(edges$from, edges$to), nodes$id)) != 0) {
    stop(
      "Nodes from edges df are missing in nodes df: ",
      paste(setdiff(union(edges$from, edges$to), nodes$id), collapse = "; ")
    )
  }
  
  gf <- igraph::graph_from_data_frame(d = edges, directed = TRUE, vertices = nodes)
  
  return(gf)
}


################################################################################

