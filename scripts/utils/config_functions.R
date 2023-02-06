


#' get sample metadata
#'
#' @param file metadata file
#'
#' @return a sample metadata tibble
#' @export
#'
#' @examples
get_metadata <- function(file){
  genus <- confs$genus
  genusPattern <- paste("(", genus, " )(?!sp\\.)", sep = "")
  
  ftype <- tools::file_ext(file)
  
  stopifnot(
    file.exists(file),
    any(ftype %in% c("csv", "tab", "tsv"))
  )
  
  if(ftype %in% c("tab", "tsv")){
    sampleInfo <- suppressMessages(readr::read_tsv(file)) 
  } else if(ftype == "csv"){
    sampleInfo <- suppressMessages(readr::read_csv(file))
  }
  
  genomeFormat <- paste("%0", stringr::str_length(max(sampleInfo$Genome)), "d", sep = "")
  
  sampleInfo %<>% dplyr::mutate(
    ## very IMP for converting Genome to a char. Int Genome will create many problems in downstream steps
    Genome = as.character(Genome),
    SpeciesName = stringi::stri_replace(
      str = SpeciesName, regex = genusPattern, replacement = "P. "
    ),
    SpeciesName = stringi::stri_replace(
      str = SpeciesName, regex = "((\\w)[^ ]+ )((\\w)[^ ]+ )(subsp\\..*)",
      replacement = "$2. $4. $5"
    ),
    nodeLabs = stringr::str_c(
      "(", stringi::stri_sprintf(format = genomeFormat, Genome), ") ",
      "(", SpeciesName,") ", sampleName,
      sep = ""
    )
  ) %>% 
    dplyr::select(sampleId, everything())
  
  return(sampleInfo)
}


################################################################################


#' Sample metadata as list
#'
#' @param df metadata data.frame
#' @param ... columns to select from metadata
#'
#' @return metadata as a named list
#' @export
#'
#' @examples
as.list_metadata <- function(df, ...){
  stopifnot(
    any(class(sampleInfo) %in% c("tbl", "data.frame")),
    nargs() > 1,
    tibble::has_name(df, "sampleId")
  )
  
  sampleInfoList <- dplyr::select(df, "sampleId", ...) %>% 
    purrr::transpose() %>% 
    purrr::set_names(nm = purrr::map(., "sampleId"))
  
  return(sampleInfoList)
}

################################################################################



