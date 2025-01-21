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


subsetLabels <- function(file = "", cluster = NULL, export.all = TRUE,
                         slice.n = "slice1", dir.out = "") {
  # Input validation
  if (!file.exists(file)) {
    stop("File does not exist.")
  }

  # Create output directory if it doesn't exist
  if (dir.out != "" && !dir.exists(dir.out)) {
    dir.create(dir.out, recursive = TRUE)
  }

  # Read data
  data <- readRDS(file)

  if (export.all) {
    # Validate cluster input
    if (is.null(cluster)) {
      stop("Cluster IDs were not specified.")
    }
    if (length(cluster) > 3) {
      stop("Number of clusters should be less than or equal to 3.")
    }
    if (!all(cluster %in% Idents(data))) {
      stop("One or more specified clusters do not exist in the data.")
    }

    # Process clusters
    clstrs <- list()
    for (i in seq_along(cluster)) {
      clstrs[[i]] <- subset(data, idents = cluster[i])

      # Ensure the slice exists
      if (!slice.n %in% names(clstrs[[i]]@images)) {
        stop(sprintf("Slice '%s' not found in the data", slice.n))
      }

      # Write coordinates
      output_file <- file.path(dir.out, paste0(cluster[i], "_coordinates.csv"))
      write.csv(clstrs[[i]]@images[[slice.n]]@coordinates,
                file = output_file,
                row.names = FALSE)
    }
    return(clstrs)
  } else {
    # Ensure the slice exists
    if (!slice.n %in% names(data@images)) {
      stop(sprintf("Slice '%s' not found in the data", slice.n))
    }

    # Write all coordinates
    output_file <- file.path(dir.out, "all_coordinates.csv")
    write.csv(data@images[[slice.n]]@coordinates,
              file = output_file,
              row.names = FALSE)

    return(invisible(NULL))
  }
}
