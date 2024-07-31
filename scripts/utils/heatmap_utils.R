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
    any(isa(phy, c("phylo", "dendrogram", "hclust", "logical")))
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
    show_row_dend = TRUE,
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
      cluster_rows = clust,
      use_raster = TRUE,
      raster_quality = 3
    )

    htList <- ht_species + ht_ani
  }

  return(htList)

}

################################################################################

#' Generate species key data for plotting
#'
#' @param genomes genomes to select
#' @param speciesInfo a dataframe with two columns: `Genome` and `SpeciesName`
#' @param type Return type: `wide`: a matrix of dimension n(Genome) x n(Species),
#'  `long`: data.frame in long format
#' @param markGenomes A named list of list with following structure:
#' `list(
#' vir = list(genomes = c("12", "2", "9"), color = "red"),
#' avir = list(genomes = c("1", "11", "4", "19"), color = "green")
#' )`
#' Default: `NULL`
#'
#' @return Either a matrix or data.frame
#' @export
#'
#' @examples NA
get_species_key_data <- function(genomes, speciesInfo, type = "wide", markGenomes = NULL) {
  stopifnot(
    type %in% c("wide", "long"),
    all(genomes %in% speciesInfo$genomeId),
    tibble::has_name(speciesInfo, "SpeciesName")
  )

  ## species name key heatmap
  speciesKey <- tibble::tibble(genomeId = genomes) %>%
    {
      if(!is.null(markGenomes)){
        dplyr::left_join(
          x = .,
          y = purrr::map(markGenomes, .f = "genomes") %>%
            tibble::enframe(name = "species", value = "genomeId") %>%
            tidyr::unnest(cols = genomeId),
          by = "genomeId"
        ) %>%
          tidyr::replace_na(replace = list(species = "1"))
      } else{
        dplyr::mutate(., species = "1")
      }
    } %>%
    dplyr::left_join(y = dplyr::select(speciesInfo, genomeId, SpeciesName), by = "genomeId")

  if (type == "wide") {
    speciesKey %<>%
      tidyr::pivot_wider(
        id_cols = genomeId, names_from = SpeciesName,
        values_from = species, values_fill = "0", names_sort = TRUE
      ) %>%
      tibble::column_to_rownames(var = "genomeId") %>%
      as.matrix()
  }

  return(speciesKey)
}

################################################################################

#' plot species name key as `ComplexHeatmap` heatmap
#'
#' @param genomes genomes to select
#' @param speciesInfo a dataframe with two columns: `Genome` and `SpeciesName`
#' @param markGenomes A named list of list with following structure:
#' `list(
#' vir = list(genomes = c("12", "2", "9"), color = "red"),
#' avir = list(genomes = c("1", "11", "4", "19"), color = "green")
#' )`
#' Default: `NULL`
#'
#' @return A heatmap
#' @export
#'
#' @examples NA
species_key_heatmap <- function(genomes, speciesInfo, markGenomes = NULL, ...){

  stopifnot(
    all(genomes %in% speciesInfo$genomeId),
    tibble::has_name(speciesInfo, "SpeciesName")
  )

  ## species name key heatmap
  speciesMat <- get_species_key_data(
    genomes = genomes, speciesInfo = speciesInfo,
    type = "wide", markGenomes = markGenomes
  )

  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == genomes))

  spColor <- purrr::map(markGenomes, "color") %>%
    unlist() %>%
    append(c("1" = "black", "0" = "white"))

  ht1 <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = spColor,
    # col = c("1" = "black"),
    cluster_columns = FALSE,
    column_split = 1:ncol(speciesMat), cluster_column_slices = FALSE,
    border = TRUE, column_gap = unit(0, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_title = "Species key",
    ...
  )

  return(ht1)
}

################################################################################

#' plot species name key as ggplot heatmap
#'
#' @param keyDf a data.frame
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples NA
species_key_ggplot <- function(keyDf) {
  stopifnot(
    tibble::has_name(keyDf, "SpeciesName"),
    tibble::has_name(keyDf, "Genome")
  )

  pt <- ggplot2::ggplot(
    data = keyDf,
    mapping = aes(x = SpeciesName, y = Genome), color = "black", fill = "black"
  ) +
    geom_tile() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
      plot.title = element_text(hjust = 0.5, vjust = 0)
    )

  return(pt)
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
    is.matrix(mat),
    is.null(speciesInfo) | (all(rownames(mat) %in% speciesInfo$genomeId)),
    is.null(hgAn) | all(colnames(mat) == hgAn$hg_id),
    any(isa(phy, c("phylo", "dendrogram", "hclust", "logical"))),
    is.null(markGenomes) |
      (is.list(markGenomes) & all(unlist(markGenomes) %in% rownames(mat)))
  )

  if (isa(phy, "phylo")) {
    clust <- as.dendrogram_ordered.phylo(phy = phy, sourceOrder = rownames(mat))
  } else if (isa(phy, "hclust")) {
    clust <- as.dendrogram(phy)
  } else {
    clust <- phy
  }

  if(isa(clust, "dendrogram")){
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    # https://github.com/jokergoo/ComplexHeatmap/issues/949
    stopifnot(
      # all(rownames(mat) == labels(clust)),
      setequal(rownames(mat), labels(clust))
    )
  }


  ## homology groups heatmap arguments
  ht_args <- list(
    matrix = mat,
    col = structure(
      viridisLite::viridis(n = max(3, min(max(mat), 6)) + 1, option = "B"),
      names = seq(0, max(3, min(max(mat), 7)))
    ),
    cluster_rows = clust, row_dend_reorder = FALSE,
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
