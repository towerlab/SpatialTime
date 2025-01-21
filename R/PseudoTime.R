#' PseudoTime
#' @param file Seurat object with module annotation
#' @param assay Seurat Assay
#' @param min_expr Minimum gene expression
#' @param min_cells minimum cells expression
#' @param mean_expr Mean gene expression
#' @param pvalue pvalue threashold to filter genes out
#' @param cores number of CPU cores to use
#' @param heatmap_plot plot monocle heatmap
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @import monocle
#' @export


PseudoTime <- function(file = NULL, assay = "RNA", min_expr = 0.1, min_cells = 3,
                       mean_expr = 0.1, pvalue = 0.05, cores = 4, return_obj = F) {

  if (!is(file, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  if (!assay %in% names(file@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  data <- as(as.matrix(file[[assay]]@data), 'sparseMatrix')
  pd <- new('AnnotatedDataFrame', data = file@meta.data)
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  HSMM <- newCellDataSet(data, phenoData = pd, featureData = fd, lowerDetectionLimit = 0.5, expressionFamily = negbinomial.size())

  HSMM <- estimateSizeFactors(HSMM)
  HSMM <- estimateDispersions(HSMM)

  HSMM <- detectGenes(HSMM, min_expr = min_expr)
  expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= min_cells))

  disp_table <- dispersionTable(HSMM)
  ordering_genes <- subset(disp_table, mean_expression >= mean_expr & 2*dispersion_empirical >= dispersion_fit)$gene_id
  HSMM <- setOrderingFilter(HSMM, ordering_genes)

  HSMM <- reduceDimension(HSMM, max_components=2, auto_param_selection = F)
  HSMM <- orderCells(HSMM, reverse=FALSE)

  HSMM@phenoData@data[["Pseudotime"]]=HSMM@phenoData@data[["st"]]
  HSMM_expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= min_cells))
  HSMM_filtered <- HSMM[HSMM_expressed_genes,]
  diff_test_res <- differentialGeneTest(HSMM_filtered, fullModelFormulaStr="~sm.ns(st)", cores = cores)

  diff=diff_test_res[,c("gene_short_name", "pval", "qval","use_for_ordering")]
  diff=diff[diff$pval < pvalue, ]
  diff=diff[order(diff$qval),]

  genelist <- row.names(diff)

  hsmm_sub <- HSMM[genelist,]

  if (return_obj) {
    return(HSMM)

    else {
      return(hsmm_sub)
    }
  }

  #return(hsmm_sub)

}
