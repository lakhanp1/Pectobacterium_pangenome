################################################################################
## visualize ANI data: all tree
#' Plot genome comparison values (ANI/Mash) score heatmap with provided clustering
#'
#' @param mat ANI score matrix
#' @param phy A phylo object
#' @param speciesInfo Optionally, a metadata for genomes. Must have SpeciesName and
#' Genome column. In addition to the genome similarity heatmap, a species key is
#' shown as heatmap in the pangenome
#' @param ... Other arguments to `ComplexHeatmap::Heatmap()`
#'
#' @return ComplexHeatmap with layout: species key heatmap | ANI heatmap
#' @export
#'
#' @examples
plot_species_ANI_heatmap <- function(mat, phy, speciesInfo = NULL, markGenomes = NULL, ...){
  
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), colnames(mat)),
    any(isa(phy, c("phylo", "dendrogram", "hclust")))
  )
  
  if (isa(phy, "phylo")) {
    clust <- as.dendrogram_ordered.phylo(phy = phy, sourceOrder = rownames(mat))
  } else if (isa(phy, "hclust")) {
    clust <- as.dendrogram(phy)
  } else {
    clust <- phy
  }
  
  stopifnot(
    # all(rownames(mat) == labels(clust)),
    setequal(rownames(mat), labels(clust))
  )
  
  ht_args <- list(
    matrix = mat,
    cluster_rows = clust,
    cluster_columns = clust,
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE, show_column_names = FALSE,
    ...
  )
  
  ht_args$column_title <- ifelse(
    is.null(ht_args$column_title), "similarity scores", ht_args$column_title
  )
  
  ## ANI heatmap
  ht_ani <- do.call(
    what = ComplexHeatmap::Heatmap,
    args = ht_args
  )
  
  htList <- ht_ani
  
  # species key heatmap
  if(!is.null(speciesInfo)){
    stopifnot(
      all(rownames(mat) %in% speciesInfo$genomeId),
      has_name(speciesInfo, "genomeId"),
      has_name(speciesInfo, "SpeciesName")
    )
    
    ht_species <- species_key_heatmap(
      genomes = rownames(mat), speciesInfo = speciesInfo,
      markGenomes = markGenomes,
      cluster_rows = clust
    )
    
    htList <- ht_species + ht_ani
  }
  
  return(htList)
  
}

################################################################################

#' Correct the order of the dendrogram created from the `phylo` object
#' 
#' This function fixes the ComplexHeatmap issue #305. 
#' https://github.com/jokergoo/ComplexHeatmap/issues/305
#'
#' @param phy a `phylo` object
#' @param sourceOrder a `rownames()` or `colnames()` vector from the source
#' matrix to which the `order.dendrogram()` need to match
#'
#' @return A dendrogram object with corrected order
#' @export
#'
#' @examples NA
as.dendrogram_ordered.phylo <- function(phy, sourceOrder){
  
  stopifnot(
    isa(phy, "phylo"),
    setequal(sourceOrder, phy$tip.label)
  )
  
  dnd <- as.dendrogram(ape::as.hclust.phylo(phy))
  
  # handle phylo->dendrogram object in ComplexHeatmap: 
  # fix the order.dendrogram() to correctly map to the matrix rows and columns
  # https://github.com/jokergoo/ComplexHeatmap/issues/305
  
  # head(order.dendrogram(dnd))
  # head(labels(dnd))
  # sourceOrder[head(order.dendrogram(dnd))]
  
  phyloOrder <- structure(1:ape::Ntip(phy), names = sourceOrder)
  newOrder <- unname(phyloOrder[labels(dnd)])
  dendextend::order.dendrogram(dnd) <- newOrder
  
  return(dnd)
}
################################################################################
