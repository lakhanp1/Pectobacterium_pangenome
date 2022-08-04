library(tidyverse)
library(here)
library(xml2)

## parse the XML file from NCBI and convert to tabular output

file_xml <- here::here("data", "reference_data", "assembly_result.xml")

ncbi <- xml2::read_xml(file_xml)
# xml2::xml_structure(ncbi)

# xml_text(ncbi)
xml_text(xml_find_all(ncbi, xpath = "//AssemblyAccession"))

xml_find_all(ncbi, xpath = "//Synonym//Genbank")

assemblyDf <- tibble::tibble(
  BioSampleAccn = xml_text(xml_find_all(ncbi, xpath = "//BioSampleAccn")),
  AssemblyAccession = xml_text(xml_find_all(ncbi, xpath = "//AssemblyAccession")),
  AssemblyName = xml_text(xml_find_all(ncbi, xpath = "//AssemblyName")),
  Genbank = xml_text(xml_find_all(ncbi, xpath = "//Synonym//Genbank")),
  RefSeq = xml_text(xml_find_all(ncbi, xpath = "//Synonym//RefSeq"))
)


readr::write_tsv(x = assemblyDf, file = here::here("data", "reference_data", "assembly_result.xml"))

