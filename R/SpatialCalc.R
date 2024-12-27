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

SpatialCalc <- function(file = "", factor = 1, colors = NULL, tissue = NULL) {

  x <- read_csv(file)

  if (!all(c("X","Y") %in% names(x))) {
    stop("Is this a fiji output coordinates file?")
  }

  if (!(length(colors) == length(tissue))) {
    stop("Number of reference lines and names should be equal!")
  }

  if (class(factor) == "character") {

    if (file.exists(factor)) {
      json_data <- fromJSON(factor)
      factor <- json_data$tissue_hires_scalef

    }
  }

  for (i in 1:length(colors)) {
    x <- x %>%
      arrange(desc(get(colors[i]))) %>%
      mutate(!!paste0(tissue[i], "_row") := case_when(get(colors[i]) == 255 ~ Y/factor, TRUE ~ 0),
             !!paste0(tissue[i], "_col") := case_when(get(colors[i]) == 255 ~ X/factor, TRUE ~ 0))
  }

  dist <- x %>%
    filter_at(vars(ends_with("_row"), ends_with("_col")), any_vars(. > 0)) %>%
    as.data.frame()

  return(dist)
}
