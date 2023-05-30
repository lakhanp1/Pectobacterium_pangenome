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
homology_groups_extract <- function(file, genomes, groups = "all", pav = FALSE){
  
  stopifnot(
    groups %in% c("core", "accessory", "unique", "all")
  )
  
  ## homology groups PAV matrix
  hgs <- suppressMessages(readr::read_csv(file)) %>% 
    dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
    dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "Genome_", "")) %>% 
    dplyr::rename_with(.fn = ~tolower(.x)) %>% 
    dplyr::mutate(homology_group_id = as.character(homology_group_id)) %>% 
    dplyr::select(hg = homology_group_id, tidyselect::all_of(genomes))
  
  ## binary matrix for homology_group PAV
  hgBinaryMat <- dplyr::mutate(
    hgs,
    dplyr::across(
      .cols = tidyselect::all_of(genomes),
      .fns = ~dplyr::if_else(.x > 1, true = 1, false = .x))
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
homology_groups_mat <- function(pandb, type, groups = NULL){
  type <- match.arg(arg = tolower(type), choices = c("pav", "cnv"))
  
  stopifnot(
    isa(pandb, c("OrgDb"))
  )
  
  hgs <- AnnotationDbi::select(
    x = pandb, keys = keys(pandb), columns = c("genome", "mRNA_key")
  )
  
  if(type == "cnv"){
    hgWide <- tidyr::pivot_wider(
      hgs,
      names_from = "GID",  values_from = "mRNA_key", 
      values_fn = ~length(.x), values_fill = 0
    )
  } else{
    hgWide <- tidyr::pivot_wider(
      hgs,
      names_from = "GID",  values_from = "mRNA_key", 
      values_fn = ~ length(.x[1]), values_fill = 0
    )
  }
  
  if(!is.null(groups)){
    hgWide %<>% dplyr::select(genome, tidyselect::all_of(groups))
  }
  
  hgMat <- tibble::column_to_rownames(hgWide, var = "genome") %>% 
    as.matrix()
  
  return(hgMat)
}

################################################################################
################################################################################
#' Plot homology group heatmap with provided clustering
#'
#' @param mat homology group matrix
#' @param phy A phylo object
#' @param metadata metadata for genomes. Must have SpeciesName and Genome column
#' @param width Vector of length 2 for heatmap widths
#' @param hgAn a data.frame for homology group annotation. It should have three
#' columns: `c(hg_id, hg_group, hg_label)`. Default: `NULL`
#' @param markGenomes A list with two elements named `list(compare = c(),
#' against = c())`. Default: `NULL`
#'
#' @return ComplexHeatmap with layout: species key heatmap | ANI heatmap
#' @export
#'
#' @examples
homology_group_heatmap <- function(mat, phy, metadata, width, hgAn = NULL, markGenomes = NULL){
  
  ## necessary checks
  stopifnot(
    setequal(rownames(mat), phy$tip.label),
    ## ensure the row order is same: this is because of a bug in ComplexHeatmap
    # https://github.com/jokergoo/ComplexHeatmap/issues/949
    all(rownames(mat) == phy$tip.label),
    all(rownames(mat) %in% metadata$Genome),
    is.null(hgAn) | all(colnames(mat) == hgAn$hg_id),
    any(class(phy) == "phylo"),
    length(width) == 2,
    all(is.numeric(width)),
    is.null(markGenomes) | 
      (is.list(markGenomes) & all(unlist(markGenomes) %in% rownames(mat)))
  )
  
  ## homology groups heatmap arguments
  ht_args <- list(
    matrix = hgMat,
    name = "hg",
    col = structure(
      viridisLite::viridis(n = min(3, max(hgMat))+1, option = "B"),
      names = seq(0, min(3, max(hgMat)))
    ),
    cluster_rows = ape::as.hclust.phylo(phy), row_dend_reorder = FALSE,
    cluster_columns = TRUE, cluster_column_slices = FALSE,
    show_column_names = TRUE, column_names_side = "bottom",
    column_names_gp = gpar(fontsize = 10),
    show_row_dend = TRUE, show_column_dend = FALSE,
    row_dend_width = unit(3, "cm"),
    show_row_names = FALSE,
    column_title = "Homology groups",
    width = unit(width[2], "cm")
  )
  
  if(!is.null(hgAn)){
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
  
  ## species name key heatmap
  speciesMat <- tibble::tibble(Genome = rownames(mat)) %>% 
    {
      if(!is.null(markGenomes)){
        dplyr::left_join(
          x = .,
          y = tibble::enframe(markGenomes, name = "species") %>% tidyr::unnest(cols = value),
          by = c("Genome" = "value")
        ) %>% 
          tidyr::replace_na(replace = list(species = "1"))
      } else{
        dplyr::mutate(., species = "1")
      }
    } %>% 
    dplyr::left_join(y = dplyr::select(metadata, Genome, SpeciesName), by = "Genome") %>% 
    tidyr::pivot_wider(
      id_cols = Genome, names_from = SpeciesName,
      values_from = species, values_fill = "0", names_sort = TRUE
    ) %>% 
    tibble::column_to_rownames(var = "Genome") %>% 
    as.matrix()
  
  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == phy$tip.label))
  
  ht_species <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = c("1" = "black", "0" = "white", "compare" = "red", "against" = "green"),
    # cluster_rows = as.hclust.phylo(phy), row_dend_reorder = FALSE,
    cluster_columns = FALSE,
    # column_order = speciesOrder,
    column_split = 1:ncol(speciesMat), cluster_column_slices = FALSE,
    border = TRUE, column_gap = unit(0, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_title = "Species key",
    width = unit(width[1], "cm")
  )
  
  # draw(ht_species)
  
  htList <- ht_species + ht_hg
  
  return(htList)
  
}

################################################################################

#' Get homology groups for a genomic region
#'
#' @param orgDb org.db pangenome object
#' @param genome genome number in pangenome
#' @param chr chromosome name
#' @param start start. default: `1` i.e. start of chromosome
#' @param end end position. default: `Inf` end of chromosome
#'
#' @return A vector of homology group identifiers
#' @export
#'
#' @examples
region_homology_groups <- function(orgDb, genome, chr, start = 1, end = Inf){
  
  stopifnot(
    !is.na(chr)
  )
  
  start <- ifelse(is.na(start), 1, start)
  end <- ifelse(is.na(end), Inf, end)
  
  df <- AnnotationDbi::select(
    x = orgDb, keys = genome, keytype = "genome",
    columns = c("GID", "chr_name", "start", "end")
    
  ) %>% 
    dplyr::mutate(
      dplyr::across(
        .cols = c(start, end), .fns = as.integer
      )
    ) %>% 
    dplyr::filter(chr_name == !!chr, start >= !!start, end <= !!end)
  
  return(df$GID)
}
################################################################################



