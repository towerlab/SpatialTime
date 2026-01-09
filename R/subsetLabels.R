#' subsetLabels
#' @param file Seurat object as Input
#' @param cluster Cluster names in Idents
#' @param export.all Export coordinates of all clusters in image
#' @param slice.n Image slice number
#' @param dir.out Ouput directory name
#'
#' @details
#' This function allows the selection of clusters of interest and exporting their coordinates
#'
#' @import tidyverse
#' @export

subsetLabels <- function(file = NULL, cluster = NULL, export.all = F, slice.n = "slice1", dir.out = "coord_out") {

  if (!is(file, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  if (dir.out != "" && !dir.exists(dir.out)) {
    dir.create(dir.out, recursive = TRUE)
  }

  if (!export.all) {

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

      write.csv(clstrs[[i]]@images[[slice.n]]@coordinates, file.path(dir.out, paste0(cluster[i], "_", names(clstrs[[i]]@images)[i], "_coordinates.csv")))
    }

    return(clstrs)

  } else {
    write.csv(file@images[[slice.n]]@coordinates,
              file.path(dir.out, "all_coordinates.csv"))
  }
}
