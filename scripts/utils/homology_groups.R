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
