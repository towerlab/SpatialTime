#'  CoordMerge
#'  @param file Seurat object as Input
#'  @param pattern Cluster names in Idents
#'
#'
#'  @details
#'  This function calculates and adds coordinates values to each line drawn in data frame.
#'
#'  @import Seurat
#'  @import tidyverse
#'  @import fs
#'  @export

#' @file - Direcotry with slices .csv file
#' @pattern - Regex pattern to select specific files within directory

CoordMerge <- function(files = "", pattern = NULL) {

  csv_files <- str_sort(fs::dir_ls(files), pattern = NULL)
  rdr <- readr::read_csv(csv_files)

  return(rdr)
}
