#' SpatialTime
#' @param file cluster coordinates file
#' @param id line as reference
#'
#'  @details
#'  This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @export

SpatialTime <- function(file = NULL, id = NULL) {

  if (!all(c("row","col","imagerow","imagecol") %in% names(file))) {
    stop("Is this a coordinate file?")
  }

  tissue_to <- file
  tissue_to <- tissue_to %>%
    mutate(barcode = X) %>%
    select(imagerow, imagecol, barcode)

  if (!is.null(id)) {
    tissue_from <- data %>%
      select(contains(id)) %>%
      filter_if(is.numeric, all_vars((.) != 0))

  } else {
    stop("id parameter is missing")
  }

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[[paste0(id, "_row")]]) ^ 2 + (imagecol - tissue_from[[paste0(id, "_col")]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}
