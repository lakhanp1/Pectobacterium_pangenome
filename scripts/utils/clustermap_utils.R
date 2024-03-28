#' Make clustermap.js input data in JSON format
#'
#' @param regions A dataframe with columns: `c("chr", "start", "end", "genomeId", "region_id")`
#' @param file Optional file name to write the clustermap data JSON string
#' @param flanking_region INT: length of flanking region around the clusters.
#' Default: 0
#' @param pandb pan.org.db object
#' @param group_colors Optional dataframe with colors for each homology group.
#'
#' @return A list object that can be written as JSON string using
#' `jsonlite::write_json()`
#' @export
#'
#' @examples NA
clustermap_data <- function(regions, file, flanking_region = 0, pandb, group_colors = NULL) {
  stopifnot(
    all(c("chr", "start", "end", "genomeId", "region_id") %in% names(regions))
  )
  
  if(!is.null(group_colors)){
    stopifnot(
      is.data.frame(group_colors),
      all(c("hg_id", "colour") %in% names(group_colors))
    )  
  }
  
  regionsLs <- split(x = regions, f = regions$region_id, drop = TRUE) %>%
    purrr::map(.f = as.list)
  
  hgStrand <- NULL
  clusterJsonDf <- NULL
  linksJsonDf <- NULL
  geneToGroup <- NULL
  slotNum <- 0
  hgFrequency <- NULL
  
  for (reg in regions$region_id) {
    regObj <- regionsLs[[reg]]
    
    if (flanking_region > 0) {
      regObj$oldStart <- regObj$start
      regObj$oldEnd <- regObj$end
      
      if (!is.na(regObj$start)) {
        regObj$start <- pmax(regObj$oldStart - flanking_region, 1)
      }
      
      if (!is.na(regObj$end)) {
        regObj$end <- regObj$end + flanking_region
      }
    }
    
    regHgs <- region_homology_groups(
      pandb = pandb, genome = regObj$genomeId,
      chr = regObj$chr, start = regObj$start, end = regObj$end,
      cols = c(
        "GID", "genePos", "chr_id", "chr_name", "start", "end", "strand",
        "mRNA_key", "genePos", "mRNA_id", "COG_description", "pfam_description"
      )
    ) %>%
      dplyr::rename(hg = GID)
    
    # summarize COG and PFAM annotations
    cog <- dplyr::select(regHgs, mRNA_key, COG_description) %>%
      dplyr::distinct() %>%
      dplyr::group_by(mRNA_key) %>%
      dplyr::summarise(COG = paste(COG_description, collapse = ";"))
    
    pfam <- dplyr::select(regHgs, mRNA_key, pfam_description) %>%
      dplyr::filter(!is.na(pfam_description)) %>%
      dplyr::distinct() %>%
      dplyr::group_by(mRNA_key) %>%
      dplyr::summarise(PFAM = paste(pfam_description, collapse = ";"))
    
    regHgs <- dplyr::select(regHgs, -COG_description, -pfam_description) %>%
      dplyr::distinct() %>%
      dplyr::left_join(y = cog, by = "mRNA_key") %>%
      dplyr::left_join(y = pfam, by = "mRNA_key") %>%
      tidyr::replace_na(replace = list(COG = "-", PFAM = "-")) %>%
      dplyr::mutate(
        strand = dplyr::if_else(strand == "-", -1, 1, 1),
        dplyr::across(.cols = c(start, end, genePos), .fns = as.integer),
        flanking = 0
      ) %>%
      dplyr::rename(
        uid = mRNA_key, label = mRNA_id, chr = chr_id
      ) %>%
      dplyr::arrange(start)
    
    # set flanking = 1 for the flanking genes
    if (flanking_region > 0) {
      flankingGenes <- which(regHgs$start < regObj$oldStart | regHgs$end > regObj$oldEnd)
      
      regHgs$flanking[flankingGenes] <- 1
    }
    
    # change the orientation of gene strand for better visualization
    thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
      as.list()
    
    if (!is.null(hgStrand)) {
      for (h in hgFrequency$hg) {
        if (!is.null(thisRegHgStrand[[h]]) & !is.null(hgStrand[[h]])) {
          if ((thisRegHgStrand[[h]] != hgStrand[[h]])) {
            regHgs <- invert_coordinates(regHgs)
            
            # update this region strands backup record after changing orientation
            thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
              as.list(thisRegHgStrand)
          }
          break
        }
      }
    } else {
      for (h in hgFrequency$hg) {
        if (!is.null(thisRegHgStrand[[h]])) {
          if (thisRegHgStrand[[h]] == -1) {
            regHgs <- invert_coordinates(regHgs)
            
            # update this region strands backup record after changing orientation
            thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
              as.list(thisRegHgStrand)
          }
          
          break
        }
      }
    }
    
    hgStrand <- thisRegHgStrand
    
    # finally, store the HG frequency information in decreasing order to
    hgFrequency <- dplyr::bind_rows(
      hgFrequency,
      dplyr::count(regHgs, hg)
    ) %>% 
      dplyr::count(hg, wt = n) %>% 
      dplyr::arrange(desc(n))
    
    # regHgs <- regHgs[1:5, ]
    
    regGenes <- dplyr::select(
      regHgs, uid, label, chr_name, chr, start, end, strand
    )
    
    # names should be a tibble type as it needs to be an {object} in JSON
    # and not a [list of {objects}]
    regGenes$names <- dplyr::select(
      regHgs, hg, COG, PFAM, genePos, genomeId
    )
    
    ## save groups information for making groups JSON
    geneToGroup <- dplyr::bind_rows(
      geneToGroup,
      dplyr::select(regHgs, hg, genes = uid, genomeId, flanking)
    )
    
    # for now there is only one loci in each cluster
    # in future, there can be multiple loci when a region is combination of
    # multiple fragments
    lociDf <- dplyr::mutate(
      regGenes,
      lociUid = chr, lociName = chr_name,
      lociStart = min(start),
      lociEnd = max(end)
    ) %>%
      tidyr::nest(
        genes = c(uid, label, names, chr, start, end, strand)
      ) %>%
      dplyr::select(
        uid = lociUid, name = lociName, start = lociStart, end = lociEnd,
        -chr_name, genes
      )
    
    clusterJsonDf <- dplyr::bind_rows(
      clusterJsonDf,
      tibble::tibble(
        uid = stringr::str_replace_all(regObj$region_id, pattern = "\\.", replacement = "_"),
        name = regObj$region_id, slot = slotNum,
        genomeId = regObj$genomeId,
        region_id = regObj$region_id,
        loci = list(lociDf)
      )
    )
    
    slotNum <- slotNum + 1
  }
  
  geneToGroup <- dplyr::mutate(
    geneToGroup,
    uid = dplyr::if_else(flanking == 1, paste(hg, "_flanking", sep = ""), hg)
  ) 
  
  ## make groups JSON
  groupsJsonDf <- dplyr::add_count(geneToGroup, hg, name = "groupFreq") %>%
    dplyr::summarise(
      genes = list(genes),
      hidden = FALSE,
      .by = c(uid, hg, flanking, groupFreq)
    ) %>%
    dplyr::arrange(desc(groupFreq)) %>%
    dplyr::select(uid, label = hg, everything()) %>%
    dplyr::arrange(desc(groupFreq))
  

  if(!is.null(group_colors)){
    # color HGs by provided colors
    groupsJsonDf <- dplyr::left_join(
      groupsJsonDf,
      dplyr::select(group_colors, hg_id, colour),
      by = c("label" = "hg_id")
    ) %>% 
      dplyr::mutate(
        colour = dplyr::if_else(
          flanking == 1, "rgb(255, 255, 255)", colour)
      )
    
  } else{
    
    # color HGs by their frequency
    hgFreqColors <- tibble::tibble(
      groupFreq = 1:max(groupsJsonDf$groupFreq),
      colour_hex = viridis::viridis(
        n = max(groupsJsonDf$groupFreq), option = "magma", direction = -1
      )
    ) %>%
      dplyr::rowwise() %>%
      dplyr::mutate(colour = paste(as.vector(col2rgb(colour_hex)), collapse = ", ")) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(colour = paste("rgb(", colour, ")", sep = "")) %>%
      dplyr::select(-colour_hex)
    
    groupsJsonDf <- dplyr::left_join(groupsJsonDf, hgFreqColors, by = "groupFreq") %>% 
      dplyr::mutate(
        colour = dplyr::if_else(flanking == 1, "rgb(255, 255, 255)", colour)
      )
  }
  
  ##### **********************
  # *** use this approach to create links as it only generate necessary links
  # *** i.e. genome 1 -> genome 2. Links reduced from ~5000 to 500.
  # *** However, there is a bug in clustermap.js which mix-up HG based coloring of genes.
  ##### **********************
  # slotLinks <- dplyr::select(clusterJsonDf, t_slot = slot, t_genome = genomeId) %>%
  #   dplyr::slice(-1) %>%
  #   dplyr::mutate(q_slot = t_slot - 1) %>%
  #   dplyr::left_join(
  #     y = dplyr::select(clusterJsonDf, q_slot = slot, q_genome = genomeId),
  #     by = "q_slot"
  #   ) %>%
  #   dplyr::left_join(
  #     y = dplyr::select(geneToGroup, hg, q_genome = genomeId, query = genes),
  #     by = "q_genome"
  #   ) %>%
  #   dplyr::left_join(
  #     y = dplyr::select(geneToGroup, hg, t_genome = genomeId, target = genes),
  #     by = c("t_genome", "hg")
  #   ) %>%
  #   dplyr::filter(if_all(.cols = c(query, target), ~ !is.na(.))) %>%
  #   tidyr::unite(col = uid, query, target, sep = "_link_", remove = F) %>%
  #   dplyr::mutate(identity = 1, similarity = 1)
  # 
  # linksJsonDf <- dplyr::select(slotLinks, uid, identity, similarity)
  # linksJsonDf$query <- dplyr::select(slotLinks, uid = query)
  # linksJsonDf$target <- dplyr::select(slotLinks, uid = target)
  
  ## make links JSON
  linksJsonDf <- purrr::pmap(
    .l = list(groupsJsonDf$uid, groupsJsonDf$genes),
    .f = function(h, g) {
      if (length(g) > 1) {
        linkCombs <- combn(x = g, m = 2) %>%
          t() %>%
          as.data.frame.matrix() %>%
          dplyr::rename(
            query = V1, target = V2
          ) %>%
          tidyr::unite(col = uid, query, target, sep = "_link_", remove = F) %>%
          dplyr::mutate(
            identity = 1,
            similarity = 1,
            hg = h
          )

        regLinks <- dplyr::select(
          linkCombs, uid, identity, similarity
        ) %>%
          dplyr::mutate(
            target = dplyr::select(linkCombs, uid = target),
            query = dplyr::select(linkCombs, uid = query)
          )

        return(regLinks)
      }
    }
  ) %>%
    purrr::list_rbind()
  
  mergedJsonList <- list(
    clusters = clusterJsonDf,
    links = linksJsonDf,
    groups = groupsJsonDf
  )
  
  if (!is.null(file)) {
    jsonlite::write_json(
      x = mergedJsonList,
      # pretty = 2,
      path = file
    )
  }
  
  return(mergedJsonList)
}

################################################################################
