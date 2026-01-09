#' SpatialCalc
#' @param file Fiji file containing drawed pixel coordinates
#' @param factor Seurat image scale factor
#' @param colors Fiji line colors
#' @param tissue Cluster identities assigned
#'
#' @details
#'Scale exported clusters coordinates to their respective drawed line colors using image factor
#'
#' @import tidyverse
#' @export


SpatialCalc <- function(file = "", factor = 1, colors = NULL, tissue = NULL) {

  x <- read_csv(file)

  if (!all(c("X","Y") %in% names(x))) {
    stop("Is this a fiji output coordinates file?")
  }

  if (is.character(factor)) {

    if (file.exists(factor)) {
      json_data <- fromJSON(factor)
      factor <- as.numeric(json_data$tissue_hires_scalef)
    }

  } else if (is.numeric(factor)) {
    factor <- as.numeric(factor)

  } else {
    stop("Factor must be either a numeric value or a path to a JSON file.")
  }

  if (!(length(colors) == length(tissue))) {
    stop("Number of reference lines and names should be equal!")
  }

  colors <- str_to_title(colors)

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
