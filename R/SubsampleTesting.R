#' Perform Testing in Subsamples
#'
#' Calculate \eqn{p}-values and log-fold changes for every gene and every cluster
#' in every subsample. Testing can be performed for global (1 vs rest) markers
#' and/or local (1 vs 1) markers. Parallelization can be achieved by setting a
#' futures plan, e.g., \code{plan(multicore)}.
#'
#' @param data A matrix- or data frame-like object containing the counts of
#'        gene-cell-pairs.
#' @param reference_identities A vector of reference cluster labels for every
#'        cell. Testing will be performed for every cluster within. Usually
#'        this is the stable clustering obtained from \emph{StabCell}.
#' @param subsampling_identities A matrix of cluster labels with rows
#'        corresponding to subsamples and columns corresponding to cells.
#'        \code{NaN} entries indicate that cells were not included in the
#'        subsample.
#' @param marker_alignment A \code{data.frame} with a columns for the subsample
#'        ID, the reference cluster ID, and the aligned subsample cluster ID.
#' @param cell_stability_score A numeric vector of cell stability scores. It
#'        will be passed to \code{test_fn()} as weights.
#' @param test_fn A function which calculates \eqn{p}-values for a subsample of
#'        the dataset. It is required to have the following correctly named
#'        arguments:
#'        \itemize{
#'          \item \code{data}: The subsample of the dataset.
#'          \item \code{identities}: A vector of cluster labels of the subsample.
#'          \item \code{weights}: The cell stability scores will be handed over
#'                for weighted testing. Can be ignored in the function body if
#'                unweighted testing is desired.
#'          \item \code{test_id1}: A single cluster label from in \code{identities}
#'                for which the test is performed.
#'          \item \code{test_id2}: A single cluster label from in \code{identities}
#'                for which the test is performed. If \code{NA} is used, a global
#'                (1 vs rest) test for cluster \code{test_id1} should be performed,
#'                otherwise a local (1 vs 1) test for clusters \code{test_id1}
#'                and \code{test_id2} should be performed.
#'          \item \code{test_args}: Additional arguments that may be used.
#'        }
#'        It should return a named vector of \eqn{p}-values with the corresponding
#'        gene as name.
#' @param test_args Additional objects that will be passed to \code{test_fn}.
#'        Defaults to \code{NULL}.
#' @param test_global A boolean flag whether to perform testing for global markers
#'        (1 vs rest).
#' @param test_local A boolean flag whether to perform testing for local markers
#'        (1 vs 1).
#' @param observations_in_columns A boolean flag indicating whether the columns
#'        of \code{data} correspond to observations (i.e. cells) or not. Defaults
#'        to \code{TRUE}.
#' @param packages A character vector of package names that are used by
#'        \code{test_fn}. Needed for parallelization with \code{doParallel}.
#'        Common choice may be \code{"Seurat"}. Defaults to \code{NULL}.
#' @param export A character vector of (custom) object names that are used by
#'        \code{test_fn} and are not covered by \code{packages}. Needed for
#'        parallelization with \code{doParallel}. Defaults to \code{NULL}.
#' @param verbose A boolean flag indication whether or not a status update should
#'        be given.Defaults to \code{TRUE}.
#'
#' @return A data frame with \eqn{p}- and log-fold change values. columns are
#'         cluster ID1, cluster ID2 (\code{NA} indicates global testing procedure),
#'         subsample, marker, p_val, and lfc. Results will be available for every
#'         cluster in \code{reference_identities}.
#'
#' @importFrom utils combn
#' @importFrom dplyr bind_cols full_join
#' @importFrom progressr progressor
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
SubsampleTesting <- function(
  data,
  reference_identities,
  subsampling_identities,
  marker_alignment,
  cell_stability_score,
  test_fn,
  test_args=NULL,
  test_global=TRUE,
  test_local=FALSE,
  observations_in_columns=TRUE,
  packages=NULL,
  export=NULL,
  verbose=TRUE
){
  if(!test_global & !test_local){ stop("At least of `test_global` and `test_local` must be true!") }

  if(!observations_in_columns){
    data <- t(data)
  }

  unique_reference_identities <- sort(unique(reference_identities))


  if (!getDoParRegistered()) registerDoSEQ()

  data_lognorm <- LogNormalizeData(data)

  # test local markers
  if(test_local){
    if(verbose) message("\ttest locally...\n")
    p <- progressor(
         nrow(subsampling_identities)
       * length(unique_reference_identities)
       * (length(unique_reference_identities) - 1)
       / 2
    )

    ij_pairs <- as.data.frame(t(combn(length(unique_reference_identities), 2)))
    colnames(ij_pairs) <- c("i", "j")
    settings <- crossing(ij_pairs, k=1:nrow(subsampling_identities))

    subsampling_p_values_local <- foreach(
      x=1:nrow(settings),
      .combine="bind_rows",
      .export=export,
      .packages=c("base", "Matrix", packages)
    ) %dorng% {
        i <- settings$i[x]  # cluster 1
        j <- settings$j[x]  # cluster 2
        k <- settings$k[x]  # subsample

        obs_idx <- !is.nan(subsampling_identities[k, ])
        subsampling_identities_k <- subsampling_identities[k, obs_idx]
        subsample <- data[ , obs_idx]
        subsample_lognorm <- data_lognorm[ , obs_idx]
        cell_stability_k <- cell_stability_score[obs_idx]


        reference_cluster_id1 <- unique_reference_identities[i]
        subsampling_cluster_id1 <- subset(
          marker_alignment,
          subsample == k & reference_cluster_id == reference_cluster_id1
        )$subsampling_cluster_id

        reference_cluster_id2 <- unique_reference_identities[j]
        subsampling_cluster_id2 <- subset(
          marker_alignment,
          subsample == k & reference_cluster_id == reference_cluster_id2
        )$subsampling_cluster_id

        p_val <- test_fn(
          data=if(observations_in_columns) subsample else t(subsample),
          identities=subsampling_identities_k,
          weights=cell_stability_k,
          test_id1=subsampling_cluster_id1,
          test_id2=subsampling_cluster_id2,
          test_args=test_args
        )
        df_p <- as.data.frame(p_val)

        lfc <- CalculateLogFoldChange(
          data=subsample_lognorm,
          identities=subsampling_identities_k,
          id1=subsampling_cluster_id1,
          id2=subsampling_cluster_id2
        )
        df_lfc <- as.data.frame(lfc)

        result <- bind_cols(
          id1=reference_cluster_id1,
          id2=reference_cluster_id2,
          subsample=k,
          full_join(
            df_p %>% mutate(marker=rownames(df_p)),
            df_lfc %>% mutate(marker=rownames(df_lfc)),
            by="marker"
          ) %>% select(marker, p_val, lfc)
        )

        p()
        return(result)
      }
  } else {
    subsampling_p_values_local <- data.frame()
  }

  # test global markers
  if(test_global){
    if(verbose) message("\ttest globally...\n")
    p <- progressor(nrow(subsampling_identities) * length(unique_reference_identities))

    settings <- crossing(
      i=1:(length(unique_reference_identities)),
      k=1:nrow(subsampling_identities)
    )

    subsampling_p_values_global <- foreach(
      x=1:nrow(settings),
      .combine="bind_rows",
      .export=export,
      .packages=c("base", "Matrix", "dplyr", packages)
    ) %dorng% {
        i <- settings$i[x]
        k <- settings$k[x]

        obs_idx <- !is.nan(subsampling_identities[k, ])
        subsampling_identities_k <- subsampling_identities[k, obs_idx]
        subsample <- data[ , obs_idx]
        subsample_lognorm <- data_lognorm[ , obs_idx]
        cell_stability_k <- cell_stability_score[obs_idx]

        reference_cluster_id1 <- unique_reference_identities[i]
        subsampling_cluster_id1 <- subset(
          marker_alignment,
          subsample == k & as.character(reference_cluster_id) == as.character(reference_cluster_id1)
        )$subsampling_cluster_id



        p_val <- test_fn(
          data=if(observations_in_columns) subsample else t(subsample),
          identities=subsampling_identities_k,
          weights=cell_stability_k,
          test_id1=subsampling_cluster_id1,
          test_id2=NULL,
          test_args=test_args
        )
        df_p <- as.data.frame(p_val)

        lfc <- CalculateLogFoldChange(
          data=subsample_lognorm,
          identities=subsampling_identities_k,
          id1=subsampling_cluster_id1,
          id2=NULL
        )
        df_lfc <- as.data.frame(lfc)

        result <- bind_cols(
          id1=reference_cluster_id1,
          id2=NA,
          subsample=k,
          full_join(
             df_p %>% mutate(marker=rownames(df_p)),
             df_lfc %>% mutate(marker=rownames(df_lfc)),
             by="marker"
          ) %>% select(marker, p_val, lfc)
        )

        p()
        return(result)
      }
  } else {
    subsampling_p_values_global <- data.frame()
  }

  # combine local & global results into one dataframe
  subsampling_p_values <- bind_rows(subsampling_p_values_global, subsampling_p_values_local)

  return(subsampling_p_values)
}



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
  data_normalized = log1p(data / col_sums * scale_factor) / log(logbase)
  return(data_normalized)
}



