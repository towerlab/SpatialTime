#' CoordMerge
#' @param files Seurat object as Input
#' @param pattern Cluster names in Idents
#'
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @import fs
#' @export

CoordMerge <- function(files = "", pattern = "") {

  csv_files <- str_sort(fs::dir_ls(files))

  selected_file <- csv_files[str_detect(csv_files, pattern)]

  rdr <- readr::read_csv(selected_file)

  if (!all(c("row","col","imagerow","imagecol") %in% names(rdr))) {
    stop("This file does not contain expected coordinates features.")
  }

  return(rdr)
}





