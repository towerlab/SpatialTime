#' SpatialVis
#' @param file Seurat Object as Input
#' @param st.calc Spatialtime values
#' @param spatial.by absolute or relative values
#' @param slice Select tissue slice
#' @param remove.na subset only tissue of spatialtime
#' @param return_obj Return object
#' SpatialVis
#' @import Seurat
#' @import tidyverse
#' @export
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.

SpatialVis <- function(file = NULL, st.calc = NULL, spatial.by = c("abs", "rel"), slice = "slice1", remove.na = F, return_obj = T) {

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

  if (remove.na) {

    sub <- subset(file, subset = st != 0)
    SpatialFeaturePlot(sub, features = "st", images = slice)

  } else {
    SpatialFeaturePlot(file, features = "st", images = slice)
  }

}


#' GeneVis
#' @param data Seurat object with pseudotime values in metadata
#' @param genes Genes to be plotted
#' @import Seurat
#' @import tidyverse
#' @export
#'
#' @details
#' Visualization of genes of interest using reference line as starting point
#'
GeneVis <- function(data = NULL, genes = NULL) {

  if (!is(data, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  if (is.null(genes)) {
    stop("Gene names must be provided!")
  }

  for (gene in genes) {
    gene_data <- FetchData(data, vars = gene)
    data[[gene]] <- gene_data

  }

  return(data)
}
