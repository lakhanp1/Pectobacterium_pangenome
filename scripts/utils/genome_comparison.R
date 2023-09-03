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
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    all(rownames(mat) == phy$tip.label),
    setequal(rownames(mat), colnames(mat)),
    any(class(phy) == "phylo")
  )
  
  
  ## ANI heatmap
  ht_ani <- ComplexHeatmap::Heatmap(
    matrix = mat,
    heatmap_legend_param = list(direction = "horizontal", legend_width = unit(3, "cm")),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    cluster_columns = ape::as.hclust.phylo(phy), column_dend_reorder = FALSE,
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE, show_column_names = FALSE,
    column_title = "similarity scores",
    ...
  )
  
  htList <- ht_ani
  
  # species key heatmap
  if(!is.null(speciesInfo)){
    stopifnot(
      all(rownames(mat) %in% speciesInfo$Genome),
      has_name(speciesInfo, "Genome"),
      has_name(speciesInfo, "SpeciesName")
    )
    
    ht_species <- species_key_heatmap(
      genomes = rownames(mat), speciesInfo = speciesInfo,
      markGenomes = markGenomes
    )
    
    htList <- ht_species + ht_ani
  }
  
  return(htList)
  
}

################################################################################
