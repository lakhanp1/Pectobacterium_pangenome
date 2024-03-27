
#' Make clustermap.js input data in JSON format
#'
#' @param regions A dataframe with columns: `c("chr", "start", "end", "genomeId", "region_id")`
#' @param file Optional file name to write the clustermap data JSON string
#' @param flanking_region INT: length of flanking region around the clusters.
#' Defaule: 0
#' @param pandb pan.org.db object 
#'
#' @return A list object that can be written as JSON string using 
#' `jsonlite::write_json()`
#' @export
#'
#' @examples NA
clustermap_data <- function(regions, file, flanking_region = 0, pandb){
  
  stopifnot(
    all(c("chr", "start", "end", "genomeId", "region_id") %in% names(regions))
  )
  regionsLs <- split(x = regions, f = regions$region_id, drop = TRUE) %>% 
    purrr::map(.f = as.list)
  
  hgStrand <- NULL
  clusterJsonDf <- NULL
  linksJsonDf <- NULL
  geneToGroup <- NULL
  slotNum <- 0
  
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
      dplyr::left_join(y = phageHgTypes, by = c("hg" = "hgId")) %>%
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
    
    # set function_category = flanking to the flanking genes added
    if (flanking_region > 0) {
      flankingGenes <- which(regHgs$start < regObj$oldStart | regHgs$end > regObj$oldEnd)
      
      regHgs$function_category[flankingGenes] <- "flanking"
      # regHgs$hg[flankingGenes] <- paste(
      #   regHgs$hg[flankingGenes], "_flanking", sep = ""
      # )
      
      regHgs$flanking[flankingGenes] <- 1
    }
    
    # change the orientation of gene strand for better visualization
    thisRegHgStrand <- dplyr::pull(regHgs, strand, name = hg) %>%
      as.list()
    
    if (!is.null(hgStrand)) {
      for (h in grpHgFreq$hgId) {
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
      for (h in grpHgFreq$hgId) {
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
    
    # regHgs <- regHgs[1:5, ]
    
    regGenes <- dplyr::select(
      regHgs, uid, label, chr_name, chr, start, end, strand
    )
    
    # names should be a tibble type as it needs to be an {object} in JSON
    # and not a [list of {objects}]
    regGenes$names <- dplyr::select(
      regHgs, COG, PFAM, function_category, hg, genePos, genomeId
    )
    
    ## save groups information for making groups JSON
    geneToGroup <- dplyr::bind_rows(
      geneToGroup,
      dplyr::select(regHgs, hg, genes = uid, flanking)
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
        region_id = regObj$region_id,
        loci = list(lociDf)
      )
    )
    
    slotNum <- slotNum + 1
  }
  
  ## make groups JSON
  groupsJsonDf <- dplyr::mutate(
    geneToGroup,
    uid = dplyr::if_else(flanking == 1, paste(hg, "_flanking", sep = ""), hg)
  ) %>%
    dplyr::add_count(hg, name = "groupFreq") %>%
    dplyr::summarise(
      genes = list(genes),
      hidden = FALSE,
      .by = c(uid, hg, flanking, groupFreq)
    ) %>%
    dplyr::arrange(desc(groupFreq)) %>%
    dplyr::select(uid, label = hg, everything()) %>%
    dplyr::left_join(y = phageHgTypes, by = c("uid" = "hgId")) %>%
    dplyr::mutate(
      function_category = dplyr::if_else(
        flanking == 1, "flanking", function_category
      )
    )
  
  # color HGs by their frequency
  hgFreqColors <- tibble::tibble(
    groupFreq = 1:max(groupsJsonDf$groupFreq),
    freq_color_hex = viridis::viridis(
      n = max(groupsJsonDf$groupFreq), option = "magma", direction = -1
    )
  ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(freq_color = paste(as.vector(col2rgb(freq_color_hex)), collapse = ", ")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(freq_color = paste("rgb(", freq_color, ")", sep = "")) %>%
    dplyr::select(-freq_color_hex)
  
  groupsJsonDf <- dplyr::left_join(
    groupsJsonDf, hgFreqColors,
    by = "groupFreq"
  ) %>%
    dplyr::left_join(hgFuncColors, by = "function_category") %>%
    dplyr::rename(colour = function_color) %>%
    dplyr::arrange(broad_function) %>%
    dplyr::arrange(desc(groupFreq))
  
  ## make links JSON
  linksJsonDf <- purrr::map2_dfr(
    .x = groupsJsonDf$uid,
    .y = groupsJsonDf$genes,
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
  )
  
  mergedJsonList <- list(
    clusters = clusterJsonDf,
    links = linksJsonDf,
    groups = groupsJsonDf
  )
  
  if(!is.null(file)){
    jsonlite::write_json(
      x = mergedJsonList,
      # pretty = 2,
      path = file
    )
  }
  
  return(mergedJsonList)
}

################################################################################
