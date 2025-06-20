#' PseudoTime
#' @param file Seurat object with module annotation
#' @param assay Seurat Assay
#' @param min_expr Minimum gene expression
#' @param min_cells minimum cells expression
#' @param mean_expr Mean gene expression
#' @param pvalue pvalue threashold to filter genes out
#' @param cores number of CPU cores to use
#' @param return_obj plot monocle heatmap
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @import monocle
#' @import Hmisc
#' @import fs
#' @export

PseudoTime <- function(file = NULL, assay = "RNA", min_expr = 0.1, min_cells = 3,
                       mean_expr = 0.1, pvalue = 0.05, cores = 1, return_obj = F) {

  if (!is(file, "Seurat")) {stop("File is not a Seurat object.")}

  if (!assay %in% names(file@assays)) {stop(paste("Assay", assay, "not found in the Seurat object."))}

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

  file@meta.data <- HSMM@phenoData@data
  genelist <- row.names(diff)
  hsmm_sub <- HSMM[genelist,]

  if (return_obj) {
    return(file)
  }

  return(hsmm_sub)

}

#' Pseudo2Time
#' @param file Seurat object with module annotation
#' @param assay Seurat Assay
#' @param min_expr Minimum gene expression
#' @param min_cells minimum cells expression
#' @param mean_expr Mean gene expression
#' @param pvalue pvalue threashold to filter genes out
#' @param cores number of CPU cores to use
#' @param return_obj plot monocle heatmap
#'
#' @details
#' This function calculates and adds coordinates values to each line drawn in data frame.
#'
#' @import Seurat
#' @import tidyverse
#' @import monocle
#' @import Hmisc
#' @import fs
#' @export

Pseudo2Time <- function(file = NULL, assay = "RNA", min_expr = 0.1, min_cells = 3,
                        mean_expr = 0.1, pvalue = 0.05, cores = 1, return_obj = F) {

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

  HSMM <- reduceDimension(HSMM, max_components=2, reduction_method = "DDRTree", auto_param_selection = F)
  HSMM <- orderCells(HSMM, reverse=FALSE)

  HSMM@phenoData@data[["Pseudotime"]]=HSMM@phenoData@data[["st"]]
  HSMM_expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= min_cells))
  HSMM_filtered <- HSMM[HSMM_expressed_genes,]
  diff_test_res <- differentialGeneTest(HSMM_filtered, fullModelFormulaStr="~sm.ns(st)", cores = cores)

  diff=diff_test_res[,c("gene_short_name", "pval", "qval","use_for_ordering")]
  diff=diff[diff$pval < pvalue, ]
  diff=diff[order(diff$qval),]

  file@meta.data <- HSMM@phenoData@data

  genelist <- row.names(diff)

  hsmm_sub <- HSMM[genelist,]

  if (return_obj) {
    return(file)
  }

  return(hsmm_sub)

}

#' PseudoM3Time
#' @param file Seurat object
#' @param assay Assays to choose
#' @param values Use pseudotime or spatialtime
#' @param q_cutoff q value cutoff
#' @param morans_cutoff morans values cutoff
#' @param cores CPU cores
#'
#' @details
#' This function calculates pseudotime using Monocle3
#'
#' @import Seurat
#' @import tidyverse
#' @import monocle3
#' @export

PseudoM3Time <- function(file = NULL, assay = c("RNA", "SCT"), values = c("pt", "st"), q_cutoff = 0.01, morans_cutoff = 0.05, cores = 4) {

  if (!is(file, "Seurat")) {
    stop("File is not a Seurat object.")
  }

  if (!assay %in% names(file@assays)) {
    stop(paste("Assay", assay, "not found in the Seurat object."))
  }

  data <- as(as.matrix(file@assays[[assay]]$data), "sparseMatrix")
  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  cds <- new_cell_data_set(data,cell_metadata = file@meta.data, gene_metadata = fData)

  cds <- preprocess_cds(cds)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  cds <- cluster_cells(cds)
  cds <- learn_graph(cds)
  cds <- order_cells(cds, reduction_method = "UMAP", verbose = F)

  cds@colData$Pseudotime <- pseudotime(cds)

  modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = cores)
  genes <- row.names(subset(modulated_genes, q_value <= q_cutoff & morans_I > morans_cutoff))

  if (values == "pt") {
    pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(cds@colData$Pseudotime)]
    pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))

  }

  if (values == "st") {
    pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(cds@colData$st)]
    pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
    pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
  } else {

    stop("values must be equal to pt or st")
  }

  return(pt.matrix)

}

#' GeneGet
#' @param var Pheatmap file from Pseudo2Time
#' @param n_clusters Number of clusters from heatmap
#' @details
#' Get genes present in each clusters identified from heatmap.
#'
#' @import Seurat
#' @import monocle
#' @import tidyverse
#' @export

GeneGet <- function(var = NULL, n_clusters = 2) {

  if (!is(var, "pheatmap")) {
    stop("File is not a pheatmap object.")
  }

  htmap <- as.data.frame(cutree(var$tree_row, k=n_clusters))
  colnames(htmap) <- "Cluster"
  htmap$Gene <- rownames(htmap)

  return(htmap)
}
