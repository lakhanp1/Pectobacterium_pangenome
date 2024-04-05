
#' Get genome ids used for a phenotype association
#'
#' @param phenotype a phenotype
#' @param confFile a phenotype association analysis config file
#'
#' @return
#' @export
#'
#' @examples
get_phenotype_association_genomes <- function(phenotype, confFile){
  
  stopifnot(
    !is.na(phenotype) | !is.null(phenotype)
  )
  
  assoConf <- suppressMessages(readr::read_tsv(confFile, col_types = "ccccc")) %>% 
    dplyr::filter(name == !!phenotype)
  
  associated <- list(
    compare = list(
      genomes = stringr::str_split_1(assoConf$compare, pattern = ","),
      color = "red"
    ),
    against = list(
      genomes = stringr::str_split_1(assoConf$against, pattern = ","),
      color = "green")
  )
  
  return(associated)
}

################################################################################

#' Get phenotype specific homology groups
#'
#' @param phenotype phenotype
#' @param panConf pangenome configuration list
#' @param save Optional file to save the phenotype specific groups.default: NULL
#'
#' @return A vector of homology group identifiers that are specific to given phenotype
#' @export
#'
#' @examples NA
phenotype_specific_groups <- function(phenotype, panConf, save = NULL){
  file_phenoAsso <- file.path(
    panConf$db$gene_classification$phenotypes$dir, phenotype, "phenotype_association.csv"
  )
  
  ## process phenotype association results
  res <- suppressMessages(
    readr::read_csv(file = file_phenoAsso, comment = "#")
  ) %>% 
    dplyr::rename_with(.fn = ~stringr::str_replace_all(.x, "( |-)", "_")) %>% 
    dplyr::rename_with(.fn = ~tolower(.x)) %>% 
    dplyr::filter(phenotype == .env$phenotype, fisher_exact_p_value != "No test")
  
  ## homology groups specific for a phenotype
  phenoSpecific <- dplyr::filter(
    res, 
    phenotype_members_absent == 0,
    other_phenotype_members_present == 0
  )
  
  phnSpeGrps <- dplyr::select(phenoSpecific, phenotype, homology_group_id) %>% 
    dplyr::group_by(phenotype) %>% 
    dplyr::summarize(
      homology_group_id = stringr::str_c(homology_group_id, collapse = ",")
    )
  
  ## save the data
  if(nrow(phenoSpecific) > 0){
    
    if(!is.null(save)){
      if(file.exists(save)){
        df <- suppressMessages(
          readr::read_tsv(save, col_names = c("phenotype", "homology_group_id"),
                          col_types = "cc")
        ) %>% 
          dplyr::filter(phenotype != .env$phenotype) %>% 
          dplyr::bind_rows(phnSpeGrps)
        
        readr::write_tsv(x = df, file = save, col_names = F)
      } else{
        readr::write_tsv(x = phnSpeGrps, file = save, col_names = F)
      }
    }
    
  }
  
  return(phenoSpecific$homology_group_id)
  # return(df)
}

################################################################################



