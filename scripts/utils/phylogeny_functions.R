

#' Add generic annotations to a ggtree
#'
#' @param pt a ggtree plot object
#' @param offset offset
#'
#' @return a ggtree plot object with annotation layers
#' @export
#'
#' @examples
annotate_ggtree <- function(pt, offset){
  pt2 <- pt +
    ggtreeExtra::geom_fruit(
      mapping = aes(starshape = type_material),
      geom = "geom_star", fill = "#62BFED", size = 2, starstroke=0.1,
      offset = offset, pwidth = 0.01
    ) +
    ggstar::scale_starshape_manual(
      values = c("type strain" = 1)
    ) + 
    ggnewscale::new_scale_color() +
    ## genome source
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "source", color = source),
      geom = "geom_point", shape = 15,
      pwidth = 0.01
    ) +
    scale_color_manual(
      values = c("NCBI" = "#406495", "WUR" = "#468e30", "NVWA" = "#ff6600")
    ) +
    ggnewscale::new_scale_fill() +
    ## taxonomy check status
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, fill = taxonomy_check_status),
      geom = "geom_star", starshape = 12, size = 2, color = alpha("white", 0),
      pwidth = 0.01
    ) +
    scale_fill_manual(
      values = c("Failed" = "red", "Inconclusive" = "blue",
                 "OK" = alpha("white", 0), "Corrected" = alpha("white", 0))
    ) +
    ## collection year
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, label = collection_year),
      geom = "geom_text", size = 3
    ) +
    ggnewscale::new_scale_color() +
    ## country
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, label = geo_loc_country, color = geo_loc_country),
      geom = "geom_text", size = 3, hjust = "left",
      pwidth = 0.03
    ) +
    scale_color_manual(
      values = c("Netherlands" = "#ff6600")
    ) +
    ggnewscale::new_scale_color() +
    ## virulence PCR result
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "vir_pcr", color = virulence_pcr),
      geom = "geom_point",
      pwidth = 0.01, offset = 0.05
    ) +
    scale_color_manual(
      values = c("positive" = "red", "negative" = "green"),
      na.value = alpha("white", 0)
    ) +
    ggnewscale::new_scale_color() +
    ## virulence phenotype
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "virulence", color = virulence),
      geom = "geom_point", shape = 17, size = 2,
      pwidth = 0.01, offset = 0.01
    ) +
    scale_color_manual(
      values = c("virulent" = "red", "avirulent" = "green"),
      na.value = alpha("white", 0)
    )
}

################################################################################
#' Title
#'
#' @param file 
#' @param metadata 
#' @param name 
#'
#' @return
#' @export
#'
#' @examples
build_annotated_tree <- function(file, metadata, name, outgroup = NULL){
  
  stopifnot(
    file.exists(file)
  )
  
  ## import tree
  rawTree <- ape::read.tree(file)
  
  ## set negative length edges => 0
  rawTree$edge.length[rawTree$edge.length < 0] <- 0
  
  ## root the tree
  if(!is.null(outgroup)){
    rawTree <- ape::root(phy = rawTree, outgroup = outgroup) %>% 
      ape::ladderize()
  }
  
  ## add metadata to tree
  treeTbl <- as_tibble(rawTree) %>%
    dplyr::full_join(y = metadata, by = c("label" = "Genome")) %>%
    treeio::as.treedata()
  
  ## draw tree and add annotation
  pt_tree <- ggtree::ggtree(tr = treeTbl) +
    labs(title = "kmer distance NJ tree")
  
  ## mark outgroup
  pt_tree2 <- mark_outgroup(pt = pt_tree, otg = outgroup, column = "sampleName")
  
  ## mark species of interest
  pt_tree3 <- pt_tree2 +
    ggtree::geom_nodelab(
      mapping = aes(label = node),
      node = "internal", size = 3, hjust = 1.3
    ) +
    ggtree::geom_tiplab(
      mapping = aes(color = SpeciesName, label = nodeLabs),
      size = 3, align = TRUE, linesize = 0.5
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    scale_color_manual(
      values = setNames(
        c("red", "#E57A44", "#088734", "#088734"),
        c(sampleInfoList[[outGroup]]$SpeciesName, "P. brasiliense",
          "P. carotovorum", "P. c. subsp. carotovorum")
      ),
      breaks = NULL,
      na.value = "black"
    ) +
    ggnewscale::new_scale_color()
  
  pt_tree4 <- annotate_ggtree(pt = pt_tree3, offset = 0.25)
  
  return(
    list(ggtree = pt_tree4, rawTree = rawTree)
  )
}
################################################################################
#' mark outgroup
#'
#' @param pt ggtree plot object
#' @param otg identifier or name of the outgroup. 
#' @param column column name that has outgroup identifier. Default: lable
#' @param color color for outgroup. Default: red
#'
#' @return ggtree plot object
#' @export
#'
#' @examples
mark_outgroup <- function(pt, otg, column = "label", color = "red"){
  
  stopifnot(
    tibble::has_name(pt$data, column),
    is.element(otg, pt$data[[column]])
  )
  
  pt2 <- pt +
    geom_point(
      mapping = aes(shape = !!sym(column), color = !!sym(column)),
      size = 4
    ) +
    scale_shape_manual(name = "outgroup", values = setNames(16, otg)) +
    scale_color_manual(name = "outgroup", values = setNames(color, otg)) +
    ggnewscale::new_scale_color() 
  
  return(pt2)
}
################################################################################
#' Prepare `pantools gene_classification` data for a clade in tree
#'
#' @param tr a ggtree object
#' @param node node label for the clade of interest
#' @param zoom optionally, a clade within which comparison to be restricted
#' @param name name of the phenotype. Default: pheno1
#' @param category value of the phenotype for one category. Default: Y
#'
#' @return a list with two elements: (pheno: dataframe for phenotype, include: genomes to include)
#' @export
#'
#' @examples
prepare_comparison_data <- function(tr, node, zoom = NULL, name = "pheno1", category = "Y"){
  stopifnot(
    is(object = tr, class2 = "treedata"),
    is.numeric(node),
    is.numeric(zoom) | is.null(zoom)
  )
  
  ## phenotype for genomes of interest
  phenoDf <- tidytree::offspring(.data = as_tibble(tr), .node = node, tiponly = TRUE) %>% 
    dplyr::select(Genome = label) %>% 
    dplyr::mutate(!!name := category)
  
  nodeGenomes <- stringr::str_c(phenoDf$Genome, collapse = ",")
  
  includeGenomes <- NA
  
  ## optionally negative phenotype
  if(!is.null(zoom)){
    parentSet <- tidytree::offspring(.data = as_tibble(tr), .node = zoom, tiponly = TRUE) %>% 
      dplyr::pull(label)
    
    includeGenomes <- stringr::str_c(parentSet, collapse = ",")
    
    negSet <- setdiff(parentSet, phenoDf$Genome)
    
  } else{
    negSet <- setdiff(na.omit(as_tibble(tr)$label), phenoDf$Genome)
  }
  
  phenoDf %<>% 
    dplyr::bind_rows(
      tibble::tibble(Genome = negSet, !!name := "N")
    )
  
  return(
    list(pheno = phenoDf, compare = nodeGenomes, 
         against = includeGenomes, name = name)
  )
}

################################################################################






