#' ShinySelection
#' @param file Seurat object
#' @param file.name Coordinates spots output file
#' @details
#' Shiny app for manual spots selection and coordinates export.
#'
#' @import Seurat
#' @import shiny
#' @import plotly
#' @import tidyverse
#' @import umap
#' @export

ShinySelection <- function(file = NULL, id = NULL, file.name = "selected_spots.csv") {

  spatial_coords <- GetTissueCoordinates(file)
  colnames(spatial_coords) <- c("x", "y")

  spatial_coords$barcode <- rownames(spatial_coords)

  metadata <- file@meta.data
  metadata$barcode <- rownames(metadata)

  plot_data <- merge(spatial_coords, metadata, by = "barcode")

  # Corrected color mapping
  unique_colors <- rainbow(length(unique(plot_data[[id]])))
  color_map <- setNames(unique_colors, as.character(unique(plot_data[[id]])))

  server <- function(input, output, session) {
    selected_indices <- reactiveVal(integer(0))
    view_state <- reactiveVal(NULL)

    output$plot <- renderPlotly({
      current_indices <- selected_indices()
      current_view <- view_state()

      p <- plot_ly(plot_data, x = ~x, y = ~y,
                   type = 'scatter',
                   mode = 'markers',
                   marker = list(
                     size = 8,
                     color = color_map[as.character(plot_data[[id]])],  # Direct color mapping
                     line = list(
                       color = 'black',
                       width = ifelse(seq_len(nrow(plot_data)) %in% current_indices, 2, 0)
                     )
                   ),
                   selectedpoints = current_indices - 1,  # 0-based indexing
                   source = "spatial_plot") %>%
        layout(
          clickmode = 'event+select',
          title = "Click spots to select/deselect",
          xaxis = list(title = "X coordinate", scaleanchor = "y", scaleratio = 1),
          yaxis = list(title = "Y coordinate"),
          dragmode = 'zoom'
        ) %>%
        config(
          scrollZoom = TRUE,
          displayModeBar = TRUE,
          modeBarButtonsToAdd = list(
            'pan2d',
            'zoom2d',
            'resetScale2d'
          )
        ) %>%
        event_register("plotly_click") %>%
        event_register("plotly_relayout")

      if (!is.null(current_view)) {
        p <- layout(p,
                    xaxis = list(range = current_view$xaxis$range),
                    yaxis = list(range = current_view$yaxis$range))
      }

      p
    })

    observeEvent(event_data("plotly_relayout", source = "spatial_plot"), {
      relayout_data <- event_data("plotly_relayout", source = "spatial_plot")

      if (!is.null(relayout_data)) {
        new_view <- list(
          xaxis = list(range = c(relayout_data[["xaxis.range[0]"]], relayout_data[["xaxis.range[1]"]])),
          yaxis = list(range = c(relayout_data[["yaxis.range[0]"]], relayout_data[["yaxis.range[1]"]]))
        )

        if (!any(is.na(unlist(new_view)))) {
          view_state(new_view)
        }
      }
    }, ignoreNULL = TRUE)

    observeEvent(event_data("plotly_click", source = "spatial_plot"), {
      click_data <- event_data("plotly_click", source = "spatial_plot")

      if (!is.null(click_data)) {
        current_indices <- selected_indices()
        spot_index <- click_data$pointNumber + 1

        if (spot_index %in% current_indices) {
          current_indices <- current_indices[current_indices != spot_index]
        } else {
          current_indices <- c(current_indices, spot_index)
        }

        selected_indices(current_indices)
      }
    }, ignoreNULL = TRUE)

    selected_spots_data <- reactive({
      indices <- selected_indices()
      if (length(indices) > 0) {
        plot_data[indices, c("barcode", "x", "y", id)]
      } else {
        data.frame()
      }
    })

    output$coordinates <- renderText({
      spots_data <- selected_spots_data()
      if (nrow(spots_data) > 0) {
        result <- "Selected spots:\n"

        for (i in 1:nrow(spots_data)) {
          x_val <- as.numeric(spots_data[i, "x"])
          y_val <- as.numeric(spots_data[i, "y"])
          barcode <- as.character(spots_data[i, "barcode"])
          cluster <- as.character(spots_data[i, id])

          result <- paste0(result,
                           "Barcode: ", barcode,
                           ", Cluster: ", cluster,
                           ", Coordinates: (",
                           format(x_val, nsmall = 2), ", ",
                           format(y_val, nsmall = 2), ")\n")
        }
        result
      } else {
        "No spots selected"
      }
    })

    output$export_button <- renderUI({
      actionButton("export", "Export Selected Spots")
    })

    output$clear_button <- renderUI({
      actionButton("clear", "Clear Selection")
    })

    observeEvent(input$export, {
      spots_data <- selected_spots_data()
      if (nrow(spots_data) > 0) {
        spots_data$x <- as.numeric(spots_data$x)
        spots_data$y <- as.numeric(spots_data$y)

        write.csv(spots_data, file.name, row.names = FALSE)
        print(paste("Exported", nrow(spots_data), "spots with barcodes to", file.name))
      } else {
        print("No spots selected to export")
      }
    })

    observeEvent(input$clear, {
      selected_indices(integer(0))
      print("Selection cleared")
    })

    output$selection_info <- renderText({
      paste("Spots selected:", length(selected_indices()))
    })
  }

  ui <- fluidPage(
    plotlyOutput("plot"),
    verbatimTextOutput("coordinates"),
    fluidRow(
      column(3, uiOutput("export_button")),
      column(3, uiOutput("clear_button"))
    ),
    textOutput("selection_info")
  )

  shinyApp(ui, server)
}


#` ShinySpots
#' @param seurat_obj Seurat object
#' @param coord_file Spots selection coordinates
#' @param slice.n Slice number
#' @param id Seurat object metadata feature name containing clusters identities
#' @details
#' This function allows coordinates extraction from specific spots barcodes previously selected in ShinySelection function
#'
#' @import Seurat
#' @import tidyverse
#' @export

ShinySpots <- function(seurat_obj = NULL, coord_file = NULL, slice.n = "slice1") {

  if (!is(seurat_obj, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  file <- read.csv(coord_file)
  labels <- file$barcode

  seurat_obj$barcode <- row.names(seurat_obj@meta.data)

  sub <- subset(seurat_obj, subset = barcode %in% labels)
  ref <- sub@images[[slice.n]]@coordinates

  return(ref)
}
