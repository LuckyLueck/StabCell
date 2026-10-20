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
#' @param cells_in_columns A boolean flag indicating whether the columns
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
    cells_in_columns=TRUE,
    packages=NULL,
    export=NULL
)
{

  if(!cells_in_columns){
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
        obs_idx <- c(obs_idx, sample((1:n_cell)[-obs_idx], 1))
      }
      subsample <- data[ , complementary_factor * obs_idx]



      # run cluster function
      identities <- cluster_fn(
        data=if(cells_in_columns) subsample else t(subsample),
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



#' Calculate Assignment Frequencies
#'
#' Given a the subsample clusterings and the alignment to an initial reference
#' clustering, the frequency of each observation being aligned to each cluster in
#' the reference clustering is calculated.
#'
#' @param subsampling_identities A matrix of subsampling cluster labels. Each row
#'        corresponds to a subsampling clustering. The number of columns must
#'        match the length of the reference identities.
#' @param clustering_alignment A \code{data.frame} with a columns for the subsample
#'        ID, the subsample cluster ID, and the aligned reference cluster ID.
#'
#' @return A \code{data.frame} with a rows corresponding to cells and columns
#'         corresponding to reference clusters. The values are the assignment
#'         frequencies.
#'
#' @importFrom dplyr %>% mutate row_number left_join group_by ungroup summarize n
#'             select
#' @importFrom tidyr pivot_longer pivot_wider drop_na
#'
#' @export
CalculateAssignmentFrequencies <- function(
    subsampling_identities,
    clustering_alignment
){
  subsampling_identities <- as.data.frame(subsampling_identities)

  assignment_freqs <- subsampling_identities %>%
    mutate(subsample=row_number()) %>%
    pivot_longer(-subsample, names_to="cell", values_to="subsampling_cluster_id") %>%
    drop_na() %>%
    left_join(
      clustering_alignment,
      by=c("subsample", "subsampling_cluster_id")
    ) %>%
    group_by(cell, reference_cluster_id) %>%
    summarize(n=n(), .groups="drop_last") %>%
    mutate(percent=n / sum(n)) %>%
    ungroup() %>%
    select(-n) %>%
    pivot_wider(names_from=reference_cluster_id, values_from=percent, values_fill=0.) %>%
    as.data.frame()

  # make sure observations are in correct order
  rownames(assignment_freqs) <- assignment_freqs$cell
  assignment_freqs <- assignment_freqs[colnames(subsampling_identities), -1, drop=FALSE]
  rownames(assignment_freqs) <- NULL

  return(assignment_freqs)
}


#' Normalize Stability Score
#'
#' Normalizes a stability score to be in the interval \eqn{[0, 1]}.
#'
#' @param x A numeric value corresponding to the stability score of one cell.
#' @param xmin,xmax A numeric value corresponding to the theoretical minimum and
#'        maximum value that \code{x} can have.
#'
#' @return A numeric value in the interval \eqn{[0, 1]}.
#'
#' @keywords internal
normalize <- function(x, xmin, xmax){ # normalizes x in [xmin, xmax] to the interval [0, 1]
  return((x - xmin) / (xmax - xmin))
}


#' Entropy-based Stability Score
#'
#' Calculate a single cell's stability score based on Shannon entropy.
#'
#' @param x A numeric vector of assignment probabilities of one cell to multiple
#'        clusters. Values may not be missing, must be positive and must add to 1.
#' @param normalize A boolean flag whether to normalize the score to be in the
#'        range \eqn{[0, 1]}. Defaults to \code{FALSE}.
#'
#' @return A numeric value. Higher values indicate more consistent assignment to
#'         one individual cluster.
#'
#' @importFrom stats na.exclude
#'
#' @keywords internal
neg_entropy <- function(x, normalize=FALSE){
  if(any(is.na(x))) stop("x must not contain NA")
  ne <- sum(na.exclude(x * log(x, base=2)))

  if(normalize){
    ne <- normalize(ne, xmin=-log(length(x), base=2), xmax=0)
  }

  return(ne)
}

#' Gini-based Stability Score
#'
#' Calculate a single cell's stability score based on Gini index.
#'
#' @param x A numeric vector of assignment probabilities of one cell to multiple
#'        clusters. Values may not be missing, must be positive and must add to 1.
#' @param normalize A boolean flag whether to normalize the score to be in the
#'        range \eqn{[0, 1]}. Defaults to \code{FALSE}.
#'
#' @return A numeric value. Higher values indicate more consistent assignment to
#'         one individual cluster.
#'
#' @keywords internal
neg_gini <- function(x, normalize=FALSE){
  if(any(is.na(x))) stop("x must not contain NA")
  ng <- sum(x^2) - 1

  if(normalize){
    ng <- normalize(ng, xmin=-(1-1/length(x)), xmax=0)
  }

  return(ng)
}

#' L2-norm-based Stability Score
#'
#' Calculate a single cell's stability score based on L2-norm.
#'
#' @param x A numeric vector of assignment probabilities of one cell to multiple
#'        clusters. Values may not be missing, must be positive and must add to 1.
#' @param normalize A boolean flag whether to normalize the score to be in the
#'        range \eqn{[0, 1]}. Defaults to \code{FALSE}.
#'
#' @return A numeric value. Higher values indicate more consistent assignment to
#'         one individual cluster.
#'
#' @keywords internal
L2score <- function(x, normalize=FALSE){
  if(any(is.na(x))) stop("x must not contain NA")
  s <- sqrt(sum(x^2))

  if(normalize){
    s <- normalize(s, xmin=1/sqrt(length(x)), xmax=1)
  }

  return(s)
}


#' Calculate Observation Stability Score
#'
#' Given assignment frequencies of observations to clusters, calculate how
#' consistently each observation is assigned to one individual cluster. Uniform
#' assignment frequencies yields a low scores, 100% assignment frequencies to
#' one individual cluster yields a high score.
#'
#' @param assignment_frequencies A data frame of assignment frequencies. Rows
#'        correspond to observations and columns to clusters.
#' @param metric A metric how to calculate the stability scores. Options are
#'        \code{"gini"}, \code{"entropy"}, and \code{"L2score"}. Defaults to
#'        \code{"gini"}.
#' @param normalize A boolean flag, whether the scores should be normalized to
#'        the interval \eqn{[0, 1]}. Defaults to \code{TRUE}.
#'
#' @return A vector of observation stability scores for each observation.
#'
#' @export
CalculateCellStabilityScore <- function(
    assignment_frequencies,
    metric=c("gini", "entropy", "L2score"),
    normalize=TRUE
  ){
  metric <- match.arg(metric)
  score_fn <- switch(
    metric,
    "gini" = neg_gini,
    "entropy" = neg_entropy,
    "L2score" = L2score
  )
  if(ncol(assignment_frequencies) == 1){
    # only 1 cluster available
    score <- rep(1, nrow(assignment_frequencies))
  } else {
    score <- apply(
      assignment_frequencies,
      1,
      function(x){ score_fn(x, normalize=normalize) }
    )
  }

  return(score)
}

#' Obtain Stable Clustering
#'
#' Find the cluster of highest assignment probability for each observation. Use
#' the first one in case of ties to be able to merge potentially spurious clusters.
#'
#' @param assignment_frequencies A data frame of assignment frequencies. Rows
#'        correspond to observations and columns to clusters.
#'
#' @return A vector of cluster labels for every observation.
#'
#' @export
FindStableClustering <- function(assignment_frequencies){
  max_idx <- max.col(assignment_frequencies, ties.method="first")
  identities <- names(assignment_frequencies)[max_idx]
  return(identities)
}


