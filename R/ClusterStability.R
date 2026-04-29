#' Calculate Assignment Probabilities
#'
#' Given a the subsampling clusterings and the alignment to an initial reference
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
#'         probabilities.
#'
#' @importFrom dplyr %>% mutate row_number left_join group_by ungroup summarize n
#'             select
#' @importFrom tidyr pivot_longer pivot_wider drop_na
#'
#' @export
CalculateAssignmentProbabilities <- function(
    subsampling_identities,
    clustering_alignment
){
  subsampling_identities <- as.data.frame(subsampling_identities)

  assignment_probs <- subsampling_identities %>%
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
  rownames(assignment_probs) <- assignment_probs$cell
  assignment_probs <- assignment_probs[colnames(subsampling_identities), -1, drop=FALSE]
  rownames(assignment_probs) <- NULL

  return(assignment_probs)
}


#' Normalize Stability Score
#'
#' Normalizes a stability score to be in the intervall \eqn{[0, 1]}.
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
#' Given assignment probabilities of observations to clusters, calculate how
#' consistently each observation is assigned to one individual cluster. Uniform
#' assignment probabilities yields a low scores, 100% assignment probability to
#' one individual cluster yields a high score.
#'
#' @param assignment_probabilities A data frame of assignment probabilities. Rows
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
    assignment_probabilities,
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
  if(ncol(assignment_probabilities) == 1){
    # only 1 cluster available
    score <- rep(1, nrow(assignment_probabilities))
  } else {
    score <- apply(
      assignment_probabilities,
      1,
      function(x){ score_fn(x, normalize=normalize) }
    )
  }

  return(score)
}

#' Obtain Subsample Clustering Mode
#'
#' Find the cluster of highest assignment probability for each observation. Use
#' the first one in case of ties to be able to merge potentially spurious clusters.
#'
#' @param assignment_probabilities A data frame of assignment probabilities. Rows
#'        correspond to observations and columns to clusters.
#'
#' @return A vector of cluster labels for every observation.
#'
#' @export
GetSubsampleClusteringMode <- function(assignment_probabilities){
  max_idx <- max.col(assignment_probabilities, ties.method="first")
  identities <- names(assignment_probabilities)[max_idx]
  return(identities)
}


