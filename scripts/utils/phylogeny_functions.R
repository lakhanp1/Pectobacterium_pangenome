

#' Add generic annotations to a ggtree
#'
#' @param pt a ggtree plot object
#' @param offset offset
#'
#' @return a ggtree plot object with annotation layers
#' @export
#'
#' @examples
annotate_ggtree <- function(pt, offset, annotations = NULL) {
  stopifnot(
    inherits(pt, "ggtree")
  )

  pt2 <- pt +
    ggtreeExtra::geom_fruit(
      mapping = aes(starshape = type_material),
      geom = "geom_star", fill = "#62BFED", size = 2, starstroke = 0.1,
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
    scale_color_viridis_d() +
    ggnewscale::new_scale_fill() +
    ## taxonomy check status
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, fill = taxonomy_check_status),
      geom = "geom_star", starshape = 12, size = 2, color = alpha("white", 0),
      pwidth = 0.01
    ) +
    scale_fill_manual(
      values = c(
        "Failed" = "red", "Inconclusive" = "blue", "OK" = alpha("white", 0),
        "corrected" = "#9DFBC1", "renamed" = "#046228"
      )
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
    ## virulence phenotype
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "virulence", color = virulence),
      geom = "geom_point", shape = 17, size = 2,
      pwidth = 0.01, offset = 0.05
    ) +
    scale_color_manual(
      values = c("virulent" = "red", "avirulent" = "green"),
      na.value = alpha("white", 0)
    ) +
    ggnewscale::new_scale_color() +
    ## virulence PCR result
    ggtreeExtra::geom_fruit(
      mapping = aes(y = id, x = "vir_pcr", color = virulence_pcr),
      geom = "geom_point",
      pwidth = 0.01, offset = 0.01
    ) +
    scale_color_manual(
      values = c("positive" = "red", "negative" = "green"),
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
build_annotated_tree <- function(file, metadata, name, outgroup = NULL) {
  stopifnot(
    file.exists(file)
  )

  ## import tree
  rawTree <- ape::read.tree(file)

  ## set negative length edges => 0
  rawTree$edge.length[rawTree$edge.length < 0] <- 0

  ## root the tree
  if (!is.null(outgroup)) {
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
        c(
          sampleInfoList[[outGroup]]$SpeciesName, "P. brasiliense",
          "P. carotovorum", "P. c. subsp. carotovorum"
        )
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
mark_outgroup <- function(pt, otg, column = "label", color = "red") {
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
#' Prepare `pantools gene_classification` data based on clades in phylogenetic tree
#'
#' @param tree a newick tree file path
#' @param node node labels for the clade of interest
#' @param against optionally, a vector of node labels for clade against which comparison
#' will be made
#' @param name name of the phenotype. Default: pheno1
#' @param category value of the phenotype for one category. Default: Y
#' @param type type of genomeset id: node or tip
#' @param excludeNode a vector of node identifiers to exclude additional genomes
#' from `against`
#' @param excludeTips a vector of tip labels to exclude additional genomes from `against`
#'
#' @return a list with following elements
#' \itemize{
#' \item pheno: a `data.frame` with for phenotype information two columns: `Genome`, `name`
#' \item compare: a comma separated list of tip lables or genome IDs to compare
#' \item against: a comma separated list of tip lables or genome IDs to compare against
#' \item name: name of the phenotype
#' }
#' @export
#'
#' @examples
clade_comparison_confs <- function(tree, node, type, against = NA, name, category = "Y",
                                   excludeNode = NA, excludeTips = NA) {
  tr <- ape::read.tree(tree)

  stopifnot(
    is.character(node),
    (all(is.element(node, tr$node.label)) | all(is.element(node, tr$tip.label))),
    is.na(against) |
      (is.character(against) & all(is.element(against, tr$node.label))),
    category != "N",
    is.element(type, c("node", "tip")),
    (is.na(excludeNode) | all(is.element(excludeNode, tr$node.label))),
    (is.na(excludeTips) | all(is.element(excludeTips, tr$tip.label)))
  )

  if (type == "node") {
    ## phenotype for genomes of interest
    ## IMP: use tree as tibble for tip labels as output when providing node label
    phenoDf <- tibble::tibble(Genome = character())

    for (nd in node) {
      phenoDf <- tidytree::offspring(.data = as_tibble(tr), .node = nd, tiponly = TRUE) %>%
        dplyr::select(Genome = label) %>%
        dplyr::bind_rows(phenoDf)
    }

    ## add phenotype
    phenoDf[[name]] <- category
  } else {
    ## when only a leaf (one genome) is being compared
    phenoDf <- tibble::tibble(
      Genome = node, !!name := category
    )
  }

  nodeGenomes <- stringr::str_c(phenoDf$Genome, collapse = ",")

  includeGenomes <- NA

  ## select the background against which comparison will be made
  if (all(!is.na(against))) {
    ## use another clade (either parent or independent clade) as background
    parentSet <- character()

    for (nd in against) {
      parentSet <- append(
        parentSet,
        tidytree::offspring(.data = as_tibble(tr), .node = nd, tiponly = TRUE) %>%
          dplyr::pull(label)
      )
    }

    negSet <- setdiff(parentSet, phenoDf$Genome)

    ## include only genomes from the groups being compared
    includeGenomes <- stringr::str_c(c(nodeGenomes, negSet), collapse = ",")
  } else {
    ## use all tips as background
    negSet <- setdiff(tr$tip.label, phenoDf$Genome)
  }

  ## additional exclusion of nodes and/or tips
  filterNeg <- character()

  if (all(!is.na(excludeNode))) {
    for (nd in excludeNode) {
      filterNeg <- append(
        filterNeg,
        tidytree::offspring(.data = as_tibble(tr), .node = nd, tiponly = TRUE) %>%
          dplyr::pull(label)
      )
    }
  }

  if (all(!is.na(excludeTips))) {
    filterNeg <- append(filterNeg, excludeTips)
  }

  if (length(filterNeg) > 0) {
    negSet <- setdiff(negSet, filterNeg)

    ## update includeGenomes
    includeGenomes <- stringr::str_c(c(nodeGenomes, negSet), collapse = ",")
  }


  phenoDf %<>%
    dplyr::bind_rows(tibble::tibble(Genome = negSet, !!name := "N"))

  return(
    list(
      pheno = phenoDf,
      compare = nodeGenomes,
      against = stringr::str_c(negSet, collapse = ","),
      includeSet = includeGenomes,
      name = name
    )
  )
}

################################################################################

#' Generate species key data for plotting
#'
#' @param genomes genomes to select
#' @param speciesInfo a dataframe with two columns: `Genome` and `SpeciesName`
#' @param type Return type: `wide`: a matrix of dimension n(Genome) x n(Species),
#'  `long`: data.frame in long format
#' @param markGenomes A list with two elements named `list(compare = c(),
#' against = c())`. Default: `NULL`
#'
#' @return Either a matrix or data.frame
#' @export
#'
#' @examples NA
get_species_key_data <- function(genomes, speciesInfo, type = "wide", markGenomes = NULL) {
  stopifnot(
    type %in% c("wide", "long"),
    all(genomes %in% speciesInfo$Genome),
    tibble::has_name(speciesInfo, "SpeciesName")
  )

  ## species name key heatmap
  speciesKey <- tibble::tibble(Genome = genomes) %>%
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
    dplyr::left_join(y = dplyr::select(speciesInfo, Genome, SpeciesName), by = "Genome")

  if (type == "wide") {
    speciesKey %<>%
      tidyr::pivot_wider(
        id_cols = Genome, names_from = SpeciesName,
        values_from = species, values_fill = "0", names_sort = TRUE
      ) %>%
      tibble::column_to_rownames(var = "Genome") %>%
      as.matrix()
  }

  return(speciesKey)
}
################################################################################

#' plot species name key as `ComplexHeatmap` heatmap
#'
#' @param genomes genomes to select
#' @param speciesInfo a dataframe with two columns: `Genome` and `SpeciesName`
#' @param markGenomes A list with two elements named `list(compare = c(),
#' against = c())`. Default: `NULL`
#' 
#' @return A heatmap
#' @export
#'
#' @examples NA
species_key_heatmap <- function(genomes, speciesInfo, markGenomes = NULL){
  
  stopifnot(
    all(genomes %in% speciesInfo$Genome),
    tibble::has_name(speciesInfo, "SpeciesName")
  )
  
  ## species name key heatmap
  speciesMat <- get_species_key_data(
    genomes = genomes, speciesInfo = speciesInfo,
    type = "wide", markGenomes = markGenomes
  )
  
  ## ensure the row order is same: this is because of a bug in ComplexHeatmap
  stopifnot(all(rownames(speciesMat) == genomes))
  
  ht1 <- ComplexHeatmap::Heatmap(
    matrix = speciesMat,
    name = "species_key",
    col = c("1" = "black", "0" = "white", "compare" = "red", "against" = "green"),
    cluster_columns = FALSE,
    column_split = 1:ncol(speciesMat), cluster_column_slices = FALSE,
    border = TRUE, column_gap = unit(0, "mm"),
    show_row_names = FALSE, show_column_names = TRUE,
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 12),
    column_title = "Species key"
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

#' Generate ANI heatmap
#'
#' @param aniDf a data.frame with three columns: `c("g1", "g2", "ANI")`
#'
#' @return a ggplot2 object
#' @export
#'
#' @examples NA
ani_heatmap <- function(aniDf) {
  pt <- ggplot2::ggplot(data = aniDf, mapping = aes(x = g1, y = g2)) +
    geom_tile(mapping = aes(fill = ANI)) +
    # scale_fill_viridis_c(name = "% identity", option = "B") +
    # scale_fill_gradientn() +
    scale_fill_stepsn(
      breaks = c(80, 85, 90, 91, 92, 92.5, 93, 93.5, 94, 95, 96, 97, 99),
      colors = viridisLite::viridis(n = 13, option = "B"),
      limits = c(80, 100)
    ) +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      plot.title = element_text(hjust = 0.5, vjust = 0)
    )

  return(pt)
}

################################################################################

#' Get nodepath for a `phylo` object
#'
#' @param phy `phylo` object
#' @param ...
#'
#' @return a tibble with `tip` and `nodepath` column
#' @export
#'
#' @examples NA
nodepath_df <- function(phy, ...) {
  stopifnot(
    isa(phy, "phylo")
  )

  np <- ape::nodepath(phy, ...)

  nlabs <- c(phy$tip.label, phy$node.label)
  names(nlabs) <- c(1:(length(phy$tip.label) + length(phy$node.label)))

  np <- purrr::map_dfr(
    .x = np,
    .f = function(x) {
      list(
        tip = nlabs[tail(x, n = 1)],
        nodepath = paste(nlabs[x], collapse = ";")
      )
    }
  )

  return(np)
}

################################################################################
