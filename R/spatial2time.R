#' Spatial2Time
#' @param file Seurat's cluster coordinates as input
#' @param file2 Seurat's reference cluster coordinates as input
#' @param id Export coordinates files
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import tidyverse
#' @export

Spatial2Time <- function(file = NULL, file2 = NULL, id = NULL) {

  columns <- c("row","col","imagerow","imagecol")

  for (i in c(file, file2)) {

    if (!all(columns %in% names(i))) {
    } else {

      stop("Is this a coordinate file?")
    }
  }

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



#' Spatial2TimeHD
#' @param file Seurat's cluster coordinates as input
#' @param file2 Seurat's reference cluster coordinates as input
#' @param id Export coordinates files
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import tidyverse
#' @export

Spatial2TimeHD <- function(file = NULL, file2 = NULL, id = NULL) {

  columns <- c("x","y")

  for (i in c(file, file2)) {

    if (!all(columns %in% names(i))) {
    } else {

      stop("Is this a coordinate file?")
    }
  }

  tissue_to <- file %>%
    mutate(barcode = ...1) %>%
    select(x, y, barcode)

  tissue_from <- file2 %>%
    mutate(barcode = ...1) %>%
    select(x, y, barcode)

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((x - tissue_from[["x"]]) ^ 2 + (y - tissue_from[["y"]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}



#' SpatialShinyTime
#' @param file Seurat's cluster coordinates as input
#' @param file2 Seurat coordinates of reference spots selected in ShinyApp
#' @param id Don't remember this
#' @details
#' This function calculates euclidean distances of reference spots selected in ShinyApp against cluster coordinates of interest
#' @import tidyverse
#' @export
#'
SpatialShinyTime <- function(file = NULL, file2 = NULL, id = NULL) {

  columns <- c("row","col","imagerow","imagecol")

  for (i in c(file, file2)) {

    if (!all(columns %in% names(i))) {
    } else {

      stop("Is this a coordinate file?")
    }
  }

  tissue_to <- file %>%
    mutate(barcode = ...1) %>%
    select(imagerow, imagecol, barcode)

  tissue_from <- file2 %>%
    mutate(barcode = rownames(file2)) %>%
    select(imagerow, imagecol, barcode)

  st_calc <- tissue_to %>%
    rowwise() %>%
    mutate(st_abs = min(sqrt((imagerow - tissue_from[["imagerow"]]) ^ 2 + (imagecol - tissue_from[["imagecol"]]) ^ 2))) %>%
    ungroup() %>%
    mutate(st_rel = (st_abs - min(st_abs)) / (max(st_abs) - min(st_abs)))

  return(st_calc)
}
