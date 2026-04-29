#' Find Subsample Clusters
#'
#' Generates subsamples of the data set and performs clustering on each subsample.
#' Parallelization can be achieved by setting a futures plan, e.g., \code{plan(multicore)}.
#'
#' @param data A matrix- or data frame-like object containing the counts of
#'        gene-cell-pairs.
#' @param cluster_fn A function which takes the following correctly named
#'        arguments:
#'        \itemize{
#'          \item \code{data}: The subsample of the dataset.
#'          \item \code{subsample_indices}: An integer vector specifying which
#'                cells are included in the subsample. This can be useful if other
#'                cell-specific information is used within cluster_fn. Values may
#'                be negative and indicate exclusion when using \code{complementary_pairs}.
#'          \item \code{cluster_args}: Additional arguments that may be used.
#'        }
#'        It must return a **numeric** or **factor** vector of cluster labels for
#' @param cluster_args Additional objects that will be passed to \code{cluster_fn}.
#'        Defaults to \code{NULL}.
#' @param subsampling_fraction A decimal number between 0 and 1 indicating which
#'        fraction of the data will be analyzed in the subsamples. Defaults to 0.5.
#' @param subsampling_iterations An integer specifying how often the subsampling
#'        procedure is performed. Defaults to 50.
#' @param complementary_pairs A boolean flag indicating whether or not to use
#'        complementary pairs subsampling. If \code{TRUE} (default), twice as
#'        many subsamples as specified by \code{subsampling_iterations} will be
#'        analyzed. Note, that if \code{subsampling_fraction} is not equal to
#'        0.5, the subsamples will be different in size: One will contain
#'        \code{subsampling_fraction} % and the other \code{1 - subsampling_fraction}
#'        of the data.
#' @param observations_in_columns A boolean flag indicating whether the columns
#'        of \code{data} correspond to observations (i.e. cells) or not. Defaults
#'        to \code{TRUE}.
#' @param packages A character vector of package names that are used by
#'        \code{cluster_fn}. Needed for parallelization worker setup. Common
#'        choice may be \code{"Seurat"}. Defaults to \code{NULL}.
#' @param export A character vector of (custom) object names that are used by
#'        \code{cluster_fn} and are not covered by \code{packages}. Needed for
#'        parallelization worker setup. Defaults to \code{NULL}.
#'
#' @return A matrix of cluster labels with rows corresponding to subsamples and
#'         columns corresponding to observations. \code{NaN} entries indicate
#'         that cells were not included in the subsample.
#'
#' @importFrom progressr progressor
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
SubsampleClustering <- function(
    data,
    cluster_fn,
    cluster_args=NULL,
    subsampling_fraction=0.5,
    subsampling_iterations=50,
    complementary_pairs=TRUE,
    observations_in_columns=TRUE,
    packages=NULL,
    export=NULL
)
{

  if(!observations_in_columns){
    data <- t(data)
  }

  if(complementary_pairs){
    complementary_factors <- c(1, -1)
  } else {
    complementary_factors <- c(1)
  }

  # run clustering in parallel
  p <- progressor(subsampling_iterations)
  if (!getDoParRegistered()) registerDoSEQ()
  subsampling_results <- foreach(
    iter=1:subsampling_iterations,
    .combine='rbind',
    .multicombine=TRUE,
    .packages=c("base", "stats", packages),
    .export=export
  ) %dorng% {


    # perform subsampling
    n_cell <- ncol(data)
    n_subsample <- floor(subsampling_fraction * n_cell)
    obs_idx <- sample.int(n_cell, n_subsample, replace=F)


    iter_result <- matrix(
      NaN,
      nrow=length(complementary_factors),
      ncol=ncol(data),
      dimnames=list(
        complementary_factor=complementary_factors,
        observation=NULL
      )
    )

    for(complementary_factor in complementary_factors){
      # account for complementary pairs if desired
      if(
        (complementary_factor == -1) &&
        (n_cell %% 2 == 1)
      ){
        obs_idx <- c(obs_idx, sample(1:n_cell[-obs_idx], 1))
      }
      subsample <- data[ , complementary_factor * obs_idx]



      # run cluster function
      identities <- cluster_fn(
        data=if(observations_in_columns) subsample else t(subsample),
        subsample_indices=complementary_factor * obs_idx,
        cluster_args=cluster_args
      )

      # save results
      # factor labels of identities will be disregarded in favor of level indicator here
      # --> no problem for further analysis
      iter_result[as.character(complementary_factor), complementary_factor * obs_idx] <- identities
    }

    p()
    iter_result
  }

  rownames(subsampling_results) <- NULL

  return(subsampling_results)
}

