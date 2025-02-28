#' SpatialTime
#' @param file file containing clusters coordinates scaled by image factor
#' @param reference cluster as reference
#' @param compare clusters spots to compare against reference
#' @details
#' Calculates euclidean distance of Fiji reference lines and cluster establishing a distance gradient
#'
#' @import Seurat
#' @import tidyverse
#' @export

SpatialTime <- function(file = NULL, reference = NULL, compare = NULL) {

  if (is.null(file)) {
    stop("Input data not found.")
  }

  if (!all(c("row","col","imagerow","imagecol") %in% names(compare))) {
    stop("Is this a coordinate file?")
  }

  tissue_to <- compare
  tissue_to <- tissue_to %>%
    mutate(barcode = X) %>%
    select(imagerow, imagecol, barcode)

  if (!is.null(reference)) {
    tissue_from <- file %>%
      select(contains(reference)) %>%
      filter_if(is.numeric, all_vars((.) != 0))

  } else {
    stop("id parameter is missing")
  }

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[[paste0(reference, "_row")]]) ^ 2 + (imagecol - tissue_from[[paste0(reference, "_col")]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}
