#' subsetLabels
#' @param file Seurat object as Input
#' @param cluster Cluster names in Idents
#' @param export.all Export coordinates files
#' @param slice.n Image slice number
#' @param dir.out directory output
#'
#' @details
#' This function allows the selection of clusters of interest and exporting their coordinates for further analysis
#'
#' @import Seurat
#' @import tidyverse
#' @export



subsetLabels <- function(file = NULL, cluster = NULL, export.all = TRUE,
                         slice.n = "slice1", dir.out = "coords_out") {

  if (!is(file, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  if (export.all) {
    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }

    if (length(cluster) > 3) {
      stop("Number of clusters should be less than or equal to 3.")
    }

    if (!all(cluster %in% Idents(file))) {
      stop("One or more specified clusters do not exist in the data.")
    }

    clstrs <- list()

    for (i in seq_along(cluster)) {
      clstrs[[i]] <- subset(file, idents = cluster[i])

      if (!dir.exists(dir.out)) {
        dir.create(dir.out)
      }

      write.csv(clstrs[[i]]@images[[slice.n]]@coordinates,
                file.path(dir.out, paste0(cluster[i], "_coordinates.csv")))
    }

    return(clstrs)

  } else {

    if (!dir.exists(dir.out)) {
      dir.create(dir.out)
    }

    write.csv(file@images[[slice.n]]@coordinates, file.path(dir.out, "all_coordinates.csv"))
  }
}
