#' Calculate the Log-Fold Change
#'
#' Calculate gene-wise log-fold change (lfc) of expression values between two
#' groups of cells.
#'
#' @param data A matrix- or data frame-like object containing the log-normalized
#'        counts of gene-cell-pairs.
#' @param identities A vector of cluster labels for every cell defining the grouping
#'        structure of the data.
#' @param id1 Group id of interest. Value must exist in \code{identities}.
#' @param id2 Group id of reference. Value must exist in \code{identities} or be
#'        \code{NULL} (default). If \code{NULL}, all other cells from \code{id1}
#'        are used.
#'
#' @return A named vector of lfc values with gene as names.
#'
#' @importFrom Matrix rowMeans
#'
#' @export
CalculateLogFoldChange <- function(data, identities, id1, id2=NULL){

  idx1 <- identities == id1
  idx2 <- if(is.null(id2)) identities != id1 else identities == id2

  mean1 <- rowMeans(data[ , idx1, drop=FALSE])
  mean2 <- rowMeans(data[ , idx2, drop=FALSE])
  lfc <- mean1 - mean2

  return(lfc)
}

#' log-Normalize Data
#'
#' Performs log-normalization for count data.
#'
#' @param data A matrix- or data frame-like object containing the counts data of
#'             gene-cell-pairs.
#' @param scale_factor Numeric value used to scale each cells count data. Makes
#'        the cells comparable to each other. Defaults to \code{10000}.
#' @param logbase The base of the logarithm to use. Defaults to \code{2}.
#'
#' @return The log-normalized data.
#'
#' @importFrom Matrix colSums
#'
#' @export
LogNormalizeData <- function(data, scale_factor=10000, logbase=2){
  col_sums = colSums(data)
  # data_normalized = log1p(data / col_sums * scale_factor) / log(logbase)
  data_normalized <- log1p(
    sweep(data, 2, col_sums, "/") * scale_factor
  ) / log(logbase)
  return(data_normalized)
}



