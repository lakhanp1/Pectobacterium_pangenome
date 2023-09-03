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
    x = pandb, keys = keys(pandb), columns = c("genome", "mRNA_key")
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
    hgWide %<>% dplyr::select(genome, tidyselect::all_of(groups))
  }
  
  hgMat <- tibble::column_to_rownames(hgWide, var = "genome") %>% 
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
    is.null(speciesInfo) | (all(rownames(mat) %in% speciesInfo$Genome)),
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
      all(rownames(mat) %in% speciesInfo$Genome),
      has_name(speciesInfo, "Genome"),
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
#' @param genome genome number in pangenome
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
    x = pandb, keys = genome, keytype = "genome",
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
      columns = c("genome", "chr", "chr_id", "chr_name", "start", "end", "strand")
    ) 
  ) %>% 
    dplyr::filter(genome == !!genome, chr_name == !!chr) %>% 
    dplyr::mutate(dplyr::across(c(start, end), as.integer)) %>% 
    dplyr::arrange(start)
  
  stopifnot(length(unique(hgInfo$chr_id)) == 1)
  
  chrInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = genome, keytype = "genome",
      columns = c("GID", "genome", "chr", "chr_id", "chr_name", "start", "end", "strand")
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
#' @param genome Find tandem match against a specific genome
#'
#' @return A list of homology groups where each element is a tandem group set
#' @export
#'
#' @examples
tandem_hg_match <- function(hgs, pandb, genome = NULL){
  hgInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = hgs,
      columns =  c("genome", "chr", "chr_id", "chr_name", "start", "end", "strand")
    ) 
  ) %>% 
    dplyr::mutate(dplyr::across(c(start, end), as.integer)) %>% 
    dplyr::arrange(start)
  
  
}


################################################################################

