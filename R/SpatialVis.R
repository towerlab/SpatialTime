#' @param file Seurat Object as Input
#' @param st.calc Spatialtime values
#' @param spatial.by absolute or relative values
#' @param slice Select tissue slice
#' @param return_obj Return object
#' SpatialVis
#' @import Seurat
#' @import tidyverse
#' @export
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.

SpatialVis <- function(file = NULL, st.calc = NULL, spatial.by = c("abs", "rel"), slice = "slice1", return_obj = F) {

  if (is.null(file) || is.null(st.calc)) {
    stop("Both 'file' and 'st.calc' must be provided.")
  }

  spatial.by <- match.arg(spatial.by)

  myBarcode <- rownames(file@meta.data)
  TissueID <- st.calc[match(myBarcode, st.calc$barcode), ]

  if (spatial.by == "abs") {
    file$st <- TissueID$st_abs
  } else if (spatial.by == "rel") {
    file$st <- TissueID$st_rel
  }

  file$st[is.na(file$st)] <- 0

  if (return_obj == T) {
    return(file)
  }

  SpatialFeaturePlot(file, features = "st", images = slice)

}
