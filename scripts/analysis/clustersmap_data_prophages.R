suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(skimr))
suppressPackageStartupMessages(library(ggdist))
suppressPackageStartupMessages(library(org.Pectobacterium.spp.pan.eg.db))
suppressPackageStartupMessages(library(jsonlite))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(dendextend))
suppressPackageStartupMessages(library(magrittr))
suppressPackageStartupMessages(library(ape))

rm(list = ls())

source("https://raw.githubusercontent.com/lakhanp1/omics_utils/main/RScripts/utils.R")
source("scripts/utils/config_functions.R")
source("scripts/utils/homology_groups.R")
source("scripts/utils/genome_scale_utils.R")
source("scripts/utils/heatmap_utils.R")
source("scripts/utils/phylogeny_functions.R")
source("scripts/utils/clustermap_utils.R")
################################################################################
set.seed(124)


confs <- prefix_config_paths(
  conf = suppressWarnings(configr::read.config(file = "project_config.yaml")),
  dir = "."
)

cluster_title <- "phage_grp_45"
outDir <- paste(confs$analysis$prophages$path, "/cluster_viz/", cluster_title, sep = "")
hg_color_categories <- confs$analysis$prophages$files$hg_broad_functions

# a vector of prophage identifiers that will be included in clustermap plot
region_cluster <- "phage_grp_45"
other_regions <- character(0)

subSample <- TRUE
cutHeight <- 0.5
addFlankingRegions <- TRUE
flankingRegion <- 5000

# ordering factor for prophages: "host" phylogeny, "hg_pav" for prophage HG PAV,
# "cluster_mash" for prophage MASH and "default" to use the provided
clusterOrder <- "host" # host, hg_pav, cluster_mash, default

# whether to keep custom regions at the bottom or consider during phylogeny
# based ordering
regions_phy_ordered <- TRUE

# regions to append as list of list with following structure
# list(r1 = list(chr, start, end, genomeId), r2 = list(chr, start, end, genomeId))
customRegions <- list(
  g_406_reg = list(chr = "NAK641_contig_10_consensus", start = 671040, end = 674984, genomeId = "g_406"),
  g_194_reg = list(chr = "NZ_CP059956.1", start = 3940785, end = 3958060, genomeId = "g_194")
)

regionClusters <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$files$clusters)
)

# # optional filter
# regionClusters %<>%
#   dplyr::filter(
#     SpeciesName == "P. brasiliense",
#     !prophage_id %in% c("g_408.vir_3", "g_403.vir_3", "g_399.vir_3")
#   )

treeMethod <- "kmer_upgma" # ani_upgma, kmer_upgma, core_snp_ml
hostPhy <- "kmer_upgma"

################################################################################
pangenome <- confs$data$pangenomes$pectobacterium.v2$name
panConf <- confs$data$pangenomes[[pangenome]]
panOrgDb <- org.Pectobacterium.spp.pan.eg.db

outPrefix <- paste(outDir, "/", cluster_title, sep = "")

if (!dir.exists(outDir)) {
  dir.create(outDir)
}

sampleInfo <- get_metadata(file = panConf$files$metadata, genus = confs$genus)

sampleInfoList <- as.list_metadata(
  df = sampleInfo, sampleId, sampleName, SpeciesName, strain, nodeLabs, genomeId
)

rawTree <- import_tree(
  file = confs$analysis$phylogeny[[treeMethod]]$files$tree_rooted,
  phylo = TRUE
)

coreTree <- import_tree(
  file = confs$analysis$phylogeny[[hostPhy]]$files$tree_rooted,
  phylo = TRUE
)

speciesOrder <- suppressMessages(
  readr::read_tsv(confs$analysis$phylogeny[[treeMethod]]$files$species_order)
)

## add species order factor levels to SpeciesName column
sampleInfo %<>% dplyr::mutate(
  SpeciesName = forcats::fct_relevel(SpeciesName, !!!rev(speciesOrder$SpeciesName))
)

################################################################################
## regions specific data
rawRegions <- suppressMessages(readr::read_tsv(confs$data$prophages$files$data)) %>%
  dplyr::select(prophage_id, chr, start, end, topology)

# read prophage HGs stored locally
regionHgs <- suppressMessages(
  readr::read_tsv(confs$analysis$prophages$preprocessing$files$raw_prophage_hg)
) %>%
  dplyr::select(prophage_id = id, hgs) %>%
  dplyr::mutate(
    hgs = stringr::str_split(hgs, ";")
  )

regionClusters <- dplyr::left_join(
  x = regionClusters, y = rawRegions, by = "prophage_id"
) %>%
  dplyr::left_join(y = regionHgs, by = "prophage_id")


regionList <- purrr::transpose(regionClusters) %>%
  purrr::set_names(nm = purrr::map(., "prophage_id"))

regionList <- append(regionList, customRegions)

clusterList <- dplyr::group_by(regionClusters, phage_grp) %>%
  dplyr::group_map(
    .f = ~ {
      list(
        phage_grp = .x$phage_grp[1],
        members = .x$prophage_id[.x$nFragments == 1],
        # save fragmented regions separately for now.
        fragmented = .x$prophage_id[.x$nFragments > 1],
        group_size = nrow(.x)
      )
    },
    .keep = TRUE
  ) %>%
  purrr::set_names(nm = purrr::map(., "phage_grp"))

mashTree <- ape::read.tree(
  file = confs$analysis$prophages$preprocessing$files$mash_hclust
)

hg_functions <- suppressMessages(readr::read_tsv(hg_color_categories)) %>%
  dplyr::select(hg_id = hgId, annotation, annotation_category)

function_types <- c(
  "other", sort(setdiff(hg_functions$annotation_category, c("other", "unknown"))),
  "unknown", "flanking"
)

hg_function_col <- purrr::map2(
  .x = function_types,
  .y = c(
    viridis::viridis(
      n = length(function_types) - 2, option = "turbo"
    ),
    "grey50", "white"
  ),
  .f = function(x, y) {
    tibble::tibble(
      annotation_category = x,
      colour_hex = y,
      colour = paste(
        "rgb(",
        paste(as.vector(col2rgb(y)), collapse = ", "),
        ")",
        sep = ""
      )
    )
  }
) %>%
  purrr::list_rbind()

hgColors <- dplyr::left_join(hg_functions, hg_function_col, by = "annotation_category") %>%
  dplyr::mutate(annotation = paste(annotation_category, ": ", annotation, sep = ""))
################################################################################
# prepare clusterjs JSON for a cluster/grp
region_set <- union(clusterList[[region_cluster]]$members, other_regions)

# if there are no members i.e. only show region of interests, then skip the next block
# if (!is.null(region_set) && !is.null(other_regions)) {
if (!is.null(region_set)) {
  # frequency for all HGs in the prophages in current cluster: used for heatmap plotting
  grpHgFreq <- regionList[region_set] %>%
    purrr::map("hgs") %>%
    unlist() %>%
    table() %>%
    tibble::enframe(name = "hg_id", value = "freq") %>%
    dplyr::mutate(freq = as.numeric(freq)) %>%
    dplyr::left_join(hgColors, by = "hg_id") %>%
    dplyr::arrange(desc(freq))

  # original MASH distance tree
  grpMashDnd <- ape::keep.tip(phy = mashTree, tip = region_set) %>%
    ape::as.hclust.phylo() %>%
    as.dendrogram() %>%
    dendextend::ladderize()

  plot(rev(grpMashDnd), horiz = TRUE)

  # HG PAV clustering
  hgPavMat <- purrr::map(
    .x = regionList[labels(grpMashDnd)],
    .f = function(x) {
      table(x$hgs) %>%
        tibble::enframe(name = "hgs", value = "count") %>%
        dplyr::mutate(count = as.numeric(count))
    }
  ) %>%
    purrr::list_rbind(names_to = "prophage_id") %>%
    tidyr::pivot_wider(
      id_cols = prophage_id, names_from = hgs,
      values_from = count, values_fill = 0
    ) %>%
    tibble::column_to_rownames(var = "prophage_id") %>%
    as.matrix()

  hgPavDnd <- hclust(d = dist(hgPavMat)) %>%
    as.dendrogram() %>%
    dendextend::ladderize()

  plot(rev(hgPavDnd), horiz = TRUE)
  # order.dendrogram(hgPavDnd) == sort(hgPavDnd, type = "nodes") %>% order.dendrogram() %>% rev()

  grpMemberData <- dplyr::left_join(
    x = tibble::tibble(region_id = region_set),
    y = dplyr::select(regionClusters, region_id = prophage_id, phage_grp, genomeId),
    by = "region_id"
  ) %>%
    dplyr::left_join(
      y = tibble::tibble(
        region_id = labels(grpMashDnd), mash_order = 1:dendextend::nleaves(grpMashDnd)
      ),
      by = "region_id"
    ) %>%
    dplyr::left_join(
      y = tibble::tibble(
        region_id = labels(hgPavDnd), pav_order = 1:dendextend::nleaves(hgPavDnd)
      ),
      by = "region_id"
    )

  # select only a subset of regions to visualize
  if (subSample) {
    pavDndCut <- dendextend::cutree(tree = hgPavDnd, h = cutHeight) %>%
      tibble::enframe(name = "region_id", value = "cut") %>%
      dplyr::add_count(cut, name = "count")

    grpSubset <- dplyr::slice_sample(pavDndCut, n = 1, by = cut)

    # ensure to include other_regions
    if (length(other_regions) > 0) {
      grpSubset <- dplyr::bind_rows(
        grpSubset,
        tibble::tibble(region_id = other_regions)
      ) %>%
        dplyr::distinct(region_id, .keep_all = TRUE)
    }

    grpMemberData <- dplyr::left_join(
      x = grpMemberData, y = pavDndCut, by = "region_id"
    ) %>%
      dplyr::left_join(
        y = tibble::tibble(region_id = grpSubset$region_id, shown = 1),
        by = "region_id"
      ) %>%
      tidyr::replace_na(replace = list(shown = 0))


    subHgPavDnd <- dendextend::prune(
      dend = hgPavDnd,
      leaves = setdiff(labels(hgPavDnd), grpSubset$region_id)
    ) %>%
      dendextend::ladderize()

    subMashDnd <- dendextend::prune(
      dend = grpMashDnd,
      leaves = setdiff(labels(grpMashDnd), grpSubset$region_id)
    ) %>%
      dendextend::ladderize()

    grpSubset <- dplyr::left_join(
      x = grpSubset,
      y = tibble::tibble(
        region_id = labels(subHgPavDnd),
        order = 1:nleaves(subHgPavDnd)
      ),
      by = "region_id"
    ) %>%
      dplyr::arrange(order)

    hgPavDnd <- dendextend::set(
      hgPavDnd,
      what = "labels_col",
      value = as.numeric(labels(hgPavDnd) %in% grpSubset$region_id) + 1
    )

    grpMashDnd <- dendextend::set(
      grpMashDnd,
      what = "labels_col",
      value = as.numeric(labels(grpMashDnd) %in% grpSubset$region_id) + 1
    )

    plot(subHgPavDnd, horiz = TRUE)

    # sort(subHgPavDnd, type = "nodes") %>% order.dendrogram() %>% rev() == order.dendrogram(subHgPavDnd)
    # sort(subHgPavDnd, type = "nodes") %>% labels() %>% rev() == labels(subHgPavDnd)
    # labels(subHgPavDnd)
    # labels(subHgPavDnd)[order.dendrogram(subHgPavDnd)]
    # sort(subHgPavDnd, type = "nodes") %>% labels() %>% rev()
    # order.dendrogram(subHgPavDnd)
    # sort(subHgPavDnd, type = "nodes") %>% order.dendrogram()

    slots <- rev(grpSubset$region_id)
    grpMemberData <- dplyr::filter(grpMemberData, shown == 1)
  } else {
    slots <- labels(grpMashDnd)

    grpMemberData <- dplyr::mutate(grpMemberData, shown = 1)
  }
} else {
  grpMemberData <- NULL
}


# add custom region genomeId and metadata
if (!is.null(customRegions) & is.list(customRegions)) {
  grpMemberData <- dplyr::bind_rows(
    x = grpMemberData,
    y = purrr::map(
      .x = customRegions,
      .f = function(x) {
        tibble::tibble(genomeId = x$genomeId, phage_grp = "custom_region", shown = 1)
      }
    ) %>%
      purrr::list_rbind(names_to = "region_id")
  )
}

# finally, add the host phylogeny order
grpMemberData %<>% dplyr::left_join(
  y = tibble::tibble(
    genomeId = ordered_tips_phylo(phy = coreTree),
    host_order = ape::Ntip.phylo(coreTree):1
  ),
  by = "genomeId"
) %>%
  dplyr::left_join(
    y = dplyr::select(sampleInfo, genomeId, SpeciesName, isolation_source, geo_loc_country),
    by = "genomeId"
  )

# optionally, do not order appended regions as per phylogeny and keep them at bottom
if (!regions_phy_ordered) {
  grpMemberData %<>%
    dplyr::mutate(
      host_order = dplyr::if_else(phage_grp == "custom_region", NA, host_order)
    )
}

if (is.null(grpMemberData$pav_order)) {
  grpMemberData$pav_order <- grpMemberData$host_order
}

if (is.null(grpMemberData$mash_order)) {
  grpMemberData$mash_order <- grpMemberData$host_order
}


# decide the order in which regions will be shown
viewRegions <- dplyr::case_when(
  clusterOrder == "host" ~
    grpMemberData$region_id[order(grpMemberData$host_order)],
  clusterOrder == "hg_pav" ~
    grpMemberData$region_id[order(grpMemberData$pav_order)],
  clusterOrder == "cluster_mash" ~
    grpMemberData$region_id[order(grpMemberData$mash_order)],
  clusterOrder == "default" ~ grpMemberData$region_id,
  TRUE ~ grpMemberData$region_id
)


regions_df <- purrr::map(
  .x = regionList[viewRegions],
  .f = function(x) {
    tibble::tibble(
      genomeId = x$genomeId, chr = x$chr, start = x$start, end = x$end
    )
  }
) %>%
  purrr::list_rbind(names_to = "region_id")


# create JSON data structure for clustermap.js
cmJson <- clustermap_data(
  regions = regions_df,
  flanking_region = ifelse(addFlankingRegions, flankingRegion, 0),
  pandb = panOrgDb,
  group_colors = hgColors,
  group_annotations = dplyr::select(hgColors, hg_id, annotation),
  file = paste(outPrefix, ".json", sep = "")
)

dplyr::select(regions_df, -genomeId) %>%
  dplyr::left_join(grpMemberData, by = "region_id") %>%
  dplyr::select(
    region_id, SpeciesName, isolation_source, geo_loc_country, genomeId,
    everything()
  ) %>%
  readr::write_tsv(
    file = paste(outPrefix, ".shown_clusters.txt", sep = "")
  )

################################################################################
# save the subtree of the core-snp tree with species names

regions_df %<>% dplyr::left_join(
  y = dplyr::select(sampleInfo, genomeId, SpeciesName), by = "genomeId"
) %>%
  dplyr::mutate(
    node_label = paste(genomeId, "|", SpeciesName)
  )

pt_core_tree <- ape::keep.tip(phy = coreTree, tip = regions_df$genomeId) %>%
  tidytree::as_tibble() %>%
  dplyr::left_join(y = regions_df, by = c("label" = "genomeId")) %>%
  tidytree::as.treedata() %>%
  ggtree::ggtree(ladderize = FALSE) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.5))) +
  ggtree::geom_tiplab(
    mapping = aes(label = node_label),
    align = TRUE, size = 4, fontface = "italic"
  ) +
  # scale_y_continuous(expand=c(0, 10)) +
  scale_y_reverse() +
  ggtree::geom_treescale()

ggplot2::ggsave(
  filename = paste(outPrefix, ".core_subtree.pdf", sep = ""),
  plot = pt_core_tree, width = 8, height = 8
)

################################################################################
# prepare homology group PAV matrix from pan.db

markGenomeColors <- list(
  prophages = list(
    genomes = purrr::map_chr(.x = regionList[region_set], .f = "genomeId") %>%
      unname(),
    color = "#D81B60"
  )
)

# highlight fragmented regions with specific color
if (length(clusterList[[region_cluster]]$fragmented) > 0) {
  markGenomeColors[["fragmented"]] <- list(
    genomes = purrr::map_chr(
      .x = regionList[
        clusterList[[region_cluster]]$fragmented
      ],
      .f = "genomeId"
    ) %>%
      unname(),
    color = "#FFC107"
  )
}

htSpecies <- species_key_heatmap(
  genomes = rawTree$tip.label, speciesInfo = sampleInfo,
  markGenomes = markGenomeColors
)

htSpecies@heatmap_param$width <- unit(12, "cm")

hgMat <- homology_groups_mat(
  pandb = panOrgDb, type = "pav", groups = grpHgFreq$hg_id
)

hgMat <- hgMat[rawTree$tip.label, ]

funcTypeColors <- dplyr::select(
  hg_function_col,
  at = annotation_category, fill = colour_hex
) %>%
  dplyr::distinct()

hgsInPan <- AnnotationDbi::select(
  x = panOrgDb, keys = grpHgFreq$hg_id,
  columns = c(
    "genome", "genomeId", "chr_id", "chr_name",
    "start", "end", "strand", "mRNA_id", "genePos"
  )
) %>%
  dplyr::mutate(
    dplyr::across(.cols = c(start, end, genome, genePos), .fns = as.numeric)
  ) %>%
  dplyr::mutate(length = end - start) %>%
  dplyr::arrange(genome, chr_id, start)

readr::write_tsv(
  hgsInPan,
  file = paste(outPrefix, ".hgs_pangenome_data.txt", sep = "")
)

hgLengths <- split(x = hgsInPan$length, f = hgsInPan$GID)
hgLengths <- hgLengths[grpHgFreq$hg_id]

# HG frequency barplot annotation
anFreq <- ComplexHeatmap::HeatmapAnnotation(
  length = ComplexHeatmap::anno_boxplot(
    x = hgLengths, which = "column"
  ),
  freq = ComplexHeatmap::anno_barplot(
    x = grpHgFreq$freq, which = "column"
  ),
  func = grpHgFreq$annotation_category,
  col = list(
    func = deframe(funcTypeColors)
  ),
  gp = gpar(col = "black"),
  annotation_height = c(4, 4, 1), height = unit(3, "cm"),
  gap = unit(1, "mm"),
  border = c(length = TRUE, freq = TRUE, func = TRUE)
)

# HG PAV heatmap using PAV data for clustering
ht <- homology_group_heatmap(
  mat = hgMat, phy = TRUE,
  # width = unit(10, "cm"),
  name = "hgs",
  top_annotation = anFreq,
  use_raster = TRUE, raster_quality = 3
)

ht@column_dend_param$cluster <- FALSE
ht@column_names_param$show <- TRUE

htList <- htSpecies + ht

pdf(
  file = paste(outPrefix, ".hg_heatmap.pdf", sep = ""),
  width = 16, height = 10
)

ComplexHeatmap::draw(
  object = htList,
  column_title = paste(
    nrow(grpHgFreq), "homology groups for prophages in cluster:", cluster_title
  ),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  main_heatmap = "hgs",
  row_dend_side = "left",
  merge_legends = TRUE,
  heatmap_legend_side = "right"
)

# use the phylogeny for clustering genomes
ComplexHeatmap::draw(
  object = htList,
  column_title = paste(
    nrow(grpHgFreq), "homology groups for prophages in cluster:", cluster_title
  ),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  main_heatmap = "hgs",
  row_dend_side = "left",
  merge_legends = TRUE,
  cluster_rows = as.dendrogram_ordered.phylo(
    phy = rawTree, sourceOrder = rawTree$tip.label
  ),
  heatmap_legend_side = "right"
)

dendextend::tanglegram(
  dend1 = rev(grpMashDnd), main_left = "MASH",
  dend2 = rev(hgPavDnd), main_right = "HG PAV",
  lwd = 2, edge.lwd = 2, lab.cex = 1, margin_inner = 5, margin_outer = 2,
  common_subtrees_color_branches = TRUE,
  highlight_distinct_edges = FALSE,
  highlight_branches_lwd = FALSE
)

if (subSample) {
  dendextend::tanglegram(
    dend1 = rev(subMashDnd), main_left = "MASH",
    dend2 = rev(subHgPavDnd), main_right = "HG PAV",
    lwd = 2, edge.lwd = 2, lab.cex = 1, margin_inner = 5, margin_outer = 2,
    common_subtrees_color_branches = TRUE,
    highlight_distinct_edges = FALSE,
    highlight_branches_lwd = FALSE
  )
}

dev.off()
