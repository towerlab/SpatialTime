#' SpatialCalc
#' @param file Fiji file
#' @param factor Scale factor
#' @param colors Line colors
#' @param tissue Tissue ID
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @export


SpatialCalc <- function(file = NULL, factor = 1, colors = NULL, tissue = NULL) {

  if (!all(c("X","Y") %in% names(x))) {
    stop("Fiji files not found!")
  }

  if (!(length(colors) == length(tissue))) {
    stop("Number of reference lines and names should be equal!")
  }

  for (i in 1:length(colors)) {
    file <- file %>%
      arrange(desc(get(colors[i]))) %>%
      mutate(!!paste0(tissue[i], "_row") := case_when(get(colors[i]) == 255 ~ Y/factor, TRUE ~ 0),
             !!paste0(tissue[i], "_col") := case_when(get(colors[i]) == 255 ~ X/factor, TRUE ~ 0))
  }

  dist <- file %>%
    filter_at(vars(ends_with("_row"), ends_with("_col")), any_vars(. > 0)) %>%
    as.data.frame()

  return(dist)
}
