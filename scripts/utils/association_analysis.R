



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
  
  return(
    list(compare = assoConf$compare, against = assoConf$against) %>% 
      purrr::map(.f = ~stringr::str_split_1(.x, pattern = ","))
  )
}







