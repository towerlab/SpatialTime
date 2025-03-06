#' subset2Labels
#' @param data Seurat object as Input
#' @param cluster Cluster names in Idents
#' @param export.all Export coordinates files
#' @param dir.out directory output
#'
#' @details
#' Export clusters coordinates and all their resptive slices
#'
#' @import Seurat
#' @import tidyverse
#' @export

subset2Labels <- function(data = NULL, cluster = NULL, export.all = F, dir.out = "coord_out") {

  if (is.null(data)) {
    stop("Data is not provided.")
  }

  if (!dir.exists(dir.out)) {
    dir.create(dir.out, recursive = TRUE)
  }

  if (!endsWith(dir.out, "/")) {
    dir.out <- paste0(dir.out, "/")
  }

  if (!export.all) {

    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }

    if (length(cluster) > 3) {
      stop("Number of clusters should be less than or equal to 3 at a time.")
    }

    if (all(cluster %in% Idents(data))) {
      clstrs <- list()

      for (i in seq_along(cluster)) {
        clstrs[[i]] <- subset(data, idents = cluster[i])
      }

      for (z in seq_along(clstrs)) {

        for (y in seq_along(names(clstrs[[z]]@images))) {
          file_name <- file.path(dir.out, paste0(cluster[z], "_", names(clstrs[[z]]@images)[y], "_coordinates.csv"))
          write.csv(clstrs[[z]]@images[[y]]@coordinates, file_name)
        }
      }

      return(clstrs)

    } else {

      stop("One or more specified clusters are not present in the data.")
    }

  } else {

    for (y in seq_along(names(data@images))) {

      file_name <- file.path(dir.out, "ALL_coordinates.csv")
      write.csv(data@images[[y]]@coordinates, file_name)
    }

    return(data)
  }
}
