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

#' Get homology groups for a genomic region
#'
#' @param pandb org.db pangenome object
#' @param genome genome identifier in pangenome
#' @param chr chromosome name
#' @param start start. default: `1` i.e. start of chromosome
#' @param end end position. default: `Inf` end of chromosome
#' @param cols Optionally, vector of columns from `pandb` object to return.
#' Default: `NULL` which returns columns `"GID", "chr_name", "start", "end", "strand"`
#' @param strand If -, homology groups are returned in reversed order
#' @param overlapping Logical: If `TRUE`, include the overlapping genes too.
#' Default: FALSE
#'
#' @return A vector of homology group identifiers
#' @export
#'
#' @examples
region_homology_groups <- function(
    pandb, genome, chr, start = 1, end = Inf, strand = "+",
    cols = NULL, overlapping = FALSE
) {
  stopifnot(
    !is.na(chr)
  )

  start <- ifelse(is.na(start), 1, start)
  end <- ifelse(is.na(end), Inf, end)

  df <- suppressMessages(AnnotationDbi::select(
    x = pandb, keys = genome, keytype = "genomeId",
    columns = union(c("GID", "chr_name", "start", "end", "strand"), cols)
  )) %>%
    dplyr::mutate(
      dplyr::across(
        .cols = c(start, end), .fns = as.integer
      )
    )

  # filter for genomic coordinates
  if (overlapping) {
    df <- dplyr::filter(
      df, chr_name == !!chr, end >= .env$start, start <= .env$end
    )
  } else{
    df <- dplyr::filter(
      df, chr_name == !!chr, start >= .env$start, end <= .env$end
    )
  }

  df <- dplyr::arrange(df, start) %>%
    tibble::as_tibble()

  if(strand == "-"){
    df <- dplyr::arrange(df, desc(start))
  }

  if(is.null(cols)){
    return(df$GID)
  } else{
    return(df)
  }
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


#' Match the HGs in the provided order on the genomes in pangenome
#'
#' @param hgs A vector of homology group ids
#' @param pandb org.db pangenome object
#' @param genomes Find tandem match against a specific genomes
#' @param gap Allowed gap during tandem match. Default = 0
#'
#' @return A list of homology groups where each element is a tandem group set
#' @export
#'
#' @examples
tandem_hg_match <- function(hgs, pandb, genomes = NULL, gap = 0) {

  hgInfo <- suppressMessages(
    AnnotationDbi::select(
      x = pandb, keys = hgs,
      columns = c("genomeId", "chr_id", "start", "end", "strand", "genePos")
    )
  )%>%
    dplyr::mutate(
      dplyr::across(.cols = c(start, end, genePos), .fns = as.integer),
    ) %>%
    dplyr::arrange(start) %>%
    dplyr::group_by(genomeId, chr_id) %>%
    dplyr::mutate(n = length(unique(GID))) %>%
    dplyr::ungroup() %>%
    dplyr::filter(n >= length(!!hgs)) %>%
    tidyr::nest(data = c(GID, genePos), .by = genomeId)

  dl <- purrr::map2(
    .x = purrr::set_names(hgInfo$genomeId), .y = hgInfo$data,
    .f = function(g, d){
      split(d$genePos, d$GID)
    }
  )

  if(!is.null(genomes)){
    genomes <- intersect(genomes, names(dl))
    stopifnot(length(genomes) > 0)

    dl <- dl[genomes]
  }

  stopifnot(length(dl) > 0)

  tandemMatches <- purrr::map(
    .x = dl,
    .f = function(chrHgs){

      currentHgPos <- -1
      tandem <- TRUE

      for (hg in hgs) {

        if (currentHgPos == -1) {
          currentHgPos <- chrHgs[[hg]]
          next
        }

        posDiff <- abs(currentHgPos - chrHgs[[hg]])

        if(all(posDiff != gap + 1)){
          tandem <- FALSE
          break
        } else{
          currentHgPos <- chrHgs[[hg]]
        }
      }

      return(tandem)
    }
  )

  return(
    purrr::discard(.x = tandemMatches, .p = isFALSE) %>% names()
  )
}



################################################################################

#' Get homology groups classification for a subset of genomes
#'
#' @param pandb org.db pangenome object
#' @param genomes Genomes for which homology groups needs to be extracted
#'
#' @return A data.frame
#' @export
#'
#' @examples
sub_pangenome_hgs <- function(pandb, genomes){

  stopifnot(!is.null(genomes))

  ## binary matrix for homology_group PAV
  hgBinaryMat <- homology_groups_mat(pandb = pandb, type = "pav")

  hgSum <- matrixStats::colSums2(
    x = hgBinaryMat, useNames = T,
    rows = which(rownames(hgBinaryMat) %in% genomes)
  ) %>%
    tibble::enframe(name = "hgId", value = "nGenomes") %>%
    dplyr::filter(nGenomes != 0) %>%
    dplyr::mutate(
      class = dplyr::case_when(
        nGenomes == !!length(genomes) ~ "core",
        nGenomes == 1 ~ "unique",
        nGenomes < !!length(genomes) & nGenomes > 1 ~ "accessory"
      )
    )
}

################################################################################

