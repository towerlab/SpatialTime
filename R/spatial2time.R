#' Spatial2Time
#' @param file Seurat object as Input
#' @param fil2 Cluster names in Idents
#' @param id Export coordinates files
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @export

Spatial2Time <- function(file = NULL, file2 = NULL, id = NULL) {

  tissue_to <- file %>%
    mutate(barcode = ...1) %>%
    select(imagerow, imagecol, barcode)

  tissue_from <- file2 %>%
    mutate(barcode = ...1) %>%
    select(imagerow, imagecol, barcode)

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[["imagerow"]]) ^ 2 + (imagecol - tissue_from[["imagecol"]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}
