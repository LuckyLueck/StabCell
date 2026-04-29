#' Calculate Selection Frequency Threshold
#'
#' Given the number of genes, the expected number of selected markes in each
#' subsample, the maximum tolerable per-family error rate (PFER), and a assumption
#' a threshold for the selection frequency is calculated.
#'
#' @param p A numeric value specifying the number of genes under consideration.
#' @param q The expected number of selected genes.
#' @param PFER_max A numeric value for the maximum tolerable per-family error
#'        rate (PFER).
#' @param assumption A character vector specifying which assumption(s) are made
#'        when calculating the threshold. Options are \code{"none"} for the
#'        Meinshausen & Bühlmann bound and \code{"unimodal"} & \code{"r_concavity"}
#'        for the two different bounds derived by Shah & Samworth. Defaults to
#'        \code{"none"}.
#' @param B The number of splits to generate pairs of complementary subsamples.
#'        Defaults to \code{NULL}, which is only permissible is \code{assumption == "none"}.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability. Defaults to
#'        \code{NULL}, which is only permissible is \code{assumption == "none"}.
#'        For other assumptions the ratio \eqn{q/p} will be used.
#'
#' @return A data frame with columns \code{n_subsample_markers}, \code{PFER_max},
#'         \code{assumption}, and the corresponding \code{selection_freq_thresh}.
#'         Logical column \code{exceeds_one} signals if the calculated selection
#'         frequency threshold is >1. In this case the value will be set to 1,
#'         but the PFER bound may not be valid any more.
#'
#' @seealso \code{\link{NoAssumptionSelectionFrequencyThreshold}},
#'          \code{\link{UnimodalSelectionFrequencyThreshold}},
#'          \code{\link{RConcavitySelectionFrequencyThreshold}}
#'
#'@export
CalculateSelectionFrequencyThreshold <- function(
  p,
  q,
  PFER_max,
  B=NULL,
  assumption="none",
  theta=NULL
){
  if(is.null(B) && any(assumption != "none")){
    warning("B is required if `assumption != 'none'`. Falling back to no assumption.")
    assumption <- "none"
  }
  if(is.null(theta)) theta <- q / p

  results <- array(
    NA,
    dim=c(length(q), length(PFER_max), length(assumption)),
    dimnames=list(q=q, PFER_max=PFER_max, assumption=assumption)
  )
  for(i in 1:length(PFER_max)){
    for(j in 1:length(assumption)){
      results[, i, j] <- switch (assumption[j],
        "none"=NoAssumptionSelectionFrequencyThreshold(p, q, PFER_max[i]),
        "unimodal"=UnimodalSelectionFrequencyThreshold(p, q, PFER_max[i], B, theta),
        "r_concavity"=RConcavitySelectionFrequencyThreshold(p, q, PFER_max[i], B, theta),
        stop("Unknown assumption. Cannot calculate selection frequency threshold.")
      )
    }
  }

  results <- as.data.frame.table(results, responseName="selection_freq_thresh") %>%
    mutate(
      q=as.numeric(as.character(q)),
      PFER_max=as.numeric(as.character(PFER_max)),
      exceeds_one=FALSE
    )

  exceeds_one <- results$selection_freq_thresh > 1.
  if(any(exceeds_one)){
    warning("Selection frequency thresholds > 1 detected. They will be set to 1. PFER bound may be violated!")
    results$selection_freq_thresh[exceeds_one] <- 1.
    results$exceeds_one <- exceeds_one
  }

  return(results)
}

#' Calculate the Selection Frequency Threshold without any further assumption
#'
#' @param p The total number of genes in the dataset.
#' @param q The expected number of selected genes.
#' @param PFER_max The maximum tolerable value for the PFER.
#'
#' @return A threshold value for the selection frequency.
#'
#' @seealso \code{\link{CalculateSelectionFrequencyThreshold}}
#'
#' @export
NoAssumptionSelectionFrequencyThreshold <- function(p, q, PFER_max){
  thresh <- (q^2 / (PFER_max * p) + 1) / 2
  return(thresh)
}

#' Calculate the Selection Frequency Threshold under Unimodality
#'
#' @param p The total number of genes in the dataset.
#' @param q The expected number of selected genes.
#' @param PFER_max The maximum tolerable value for the PFER.
#' @param B The number of splits to generate pairs of complementary subsamples.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability.
#'
#' @return A threshold value for the selection frequency. May return \code{Inf}
#'         if no threshold can be calculated
#'
#' @seealso \code{\link{CalculateSelectionFrequencyThreshold}}
#'
#' @export
UnimodalSelectionFrequencyThreshold <- function(p, q, PFER_max, B, theta){
  if(theta > (1 / sqrt(3))){
    warning("Unimodal assumption inappropriate (theta > 1/sqrt(3)). Return `Inf`.")
    return(Inf)
  }

  thresh_llim <- min(
    0.5 + theta^2,
    0.5 + 1 / (2*B) + 0.75 * theta^2
  )
  thresh_candidates <- 0.5 + (2:B)/(2*B)
  thresh_candidates <- thresh_candidates[thresh_candidates > thresh_llim]
  if(length(thresh_candidates) == 0){
    warning("Unimodal assumption inappropriate (all thresh_candidates <= lower limit). Return `Inf`.")
    return(Inf)
  }

  C <- rep(NA, length(thresh_candidates))
  idx <- (thresh_candidates <= 0.75)
  C[idx] <- 1 / (2 * (2 * thresh_candidates[idx] - 1 - 1/(2*B)))
  C[!idx] <- (4 * (1 - thresh_candidates[!idx] + 1/(2*B))) / (1 + 1/B)

  PFER_ulims <- C * theta * q
  thresh <- thresh_candidates[PFER_ulims <= PFER_max][1]

  if(is.na(thresh)) return(Inf) else return(thresh)
}

#' Calculate the Selection Frequency Threshold under r-concavity
#'
#' @param p The total number of genes in the dataset.
#' @param q The expected number of selected genes.
#' @param PFER_max The maximum tolerable value for the PFER.
#' @param B The number of splits to generate pairs of complementary subsamples.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability.
#'
#' @return A threshold value for the selection frequency.
#'
#' @seealso \code{\link{CalculateSelectionFrequencyThreshold}}, \code{\link{r.TailProbs}}
#'
#' @export
RConcavitySelectionFrequencyThreshold <- function(p, q, PFER_max, B, theta){

  thresh_candidates <- 1:(2*B)/(2*B)

  PFER_ulims <- p * pmin(
    c(rep(Inf, B), r.TailProbs(theta^2, B, -0.5)),
    r.TailProbs(theta, 2*B, -0.25)
  )
  thresh <- thresh_candidates[(PFER_ulims <= PFER_max) & (thresh_candidates > theta)][1]

  if(is.na(thresh)) return(Inf) else return(thresh)
}


#' Calculate tail probabilities under r-concavity
#'
#' This is a slightly modified version of the function published by Shah & Samworth
#' alongside their article "Variable selection with error control: Another look at
#' Stability Selection" (DOI: 10.1111/j.1467-9868.2011.01034.x). Their implementation
#' is available \href{https://www.statslab.cam.ac.uk/~rjs57/r_concave_tail.R}{here}
#' or \href{https://www.statslab.cam.ac.uk/~rds37/papers/r_concave_tail.R}{here}.
#' This function is meant for internal use in \code{\link{RConcavitySelectionFrequencyThreshold}}
#' only.
#'
#' @param eta ...
#' @param B The number of splits to generate pairs of complementary subsamples.
#' @param r ...
#'
#' @return The calculated tail probabilities.
#'
#' @seealso \code{\link{RConcavitySelectionFrequencyThreshold}}
#'
#' @importFrom stats optimize uniroot
#'
#' @keywords internal
r.TailProbs <- function(eta, B, r) {

  # TailProbs returns a vector with the tail probability for each \tau = ceil{B*2\eta}/B + 1/B,...,1
  # We return 1 for all \tau = 0, 1/B, ... , ceil{B*2\eta}/B
  MAXa <- 100000
  MINa <- 0.00001
  s <- 1/r
  etaB <- eta * B
  k_start <- (ceiling(2 * etaB) + 1)
  output <- rep(Inf, B)

  if(k_start > B) return(output)

  a_vec <- rep(MAXa,B)

  Find.a <- function(prev_a) {
    uniroot(Calc.a, lower=MINa, upper=prev_a,
            tol=.Machine$double.eps^0.75
            )$root
  }
  Calc.a <- function(a) {
    denom <- sum((a + 0:k)^s)
    num <- sum((0:k) * (a + 0:k)^s)
    num / denom - etaB
  }

  for(k in k_start:B){
    a_vec[k] <- Find.a(a_vec[k-1])
  }

  OptimInt <- function(a) {
    # NB this function makes use of several global variables
    num <- (k + 1 - etaB) * sum((a + 0:(t-1))^s)
    denom <- sum((k + 1 - (0:k)) * (a + 0:k)^s)
    1 - num / denom
  }

  prev_k <- k_start
  for(t in k_start:B) {
    cur_optim <- rep(0, B)
    cur_optim[B] <- OptimInt(a_vec[B])
    if (prev_k <= (B-1)) {
      for (k in prev_k:(B-1)){
        cur_optim[k] <- optimize(
          f=OptimInt,
          lower=a_vec[k+1],
          upper=a_vec[k],
          maximum=TRUE
        )$objective
      }
    }
    output[t] <- max(cur_optim)
    prev_k <- which.max(cur_optim)
  }
  return(output)
}




#' Find Stable Markers
#'
#' Given selection frequencies and a corresponding threshold, this function finds
#' globally (1 vs rest) or locally (1 vs 1) stable markers.
#'
#' @param id1,id2 Cluster IDs for which the markers are detected. If \code{id2}
#'        is NULL, globally stable markers are found.
#' @param selection_frequencies A data frame with a columns for the \eqn{p}-value
#'        \eqn{\lambda}, cluster ID1, cluster ID2 (\code{NA} means global marker
#'        for ID1), and columns corresponding to the gene names in the data set.
#'        Values are the selection frequencies.
#' @param PFER_max A numeric value for the maximum tolerable per-family error
#'        rate (PFER).
#' @param assumption A character vector specifying which assumption is made
#'        when calculating the threshold. Options are \code{"none"} for the
#'        Meinshausen & Bühlmann bound and \code{"unimodal"} & \code{"r_concavity"}
#'        for the two different bounds derived by Shah & Samworth. Defaults to
#'        \code{"none"}.
#' @param B The number of splits to generate pairs of complementary subsamples.
#'        Defaults to \code{NULL}, which is only permissible when \code{assumption == "none"}.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability. Defaults to
#'        \code{NULL}, which is only permissible when \code{assumption == "none"}.
#'        For other assumptions the ratio \eqn{q/p} will be used.
#'
#' @return A data frame with columns for cluster IDs, \eqn{p}-value threshold,
#'         selection frequency threshold, marker, and observed selection frequency.
#'         \code{NA} entries mean no marker was detected for this combination of
#'         cluster, \eqn{p}-value threshold and selection frequency threshold.
#'
#' @seealso \code{\link{FindAllGlobalStableMarkers}}, \code{\link{FindAllLocalStableMarkers}}
#'
#' @importFrom methods as
#' @importFrom dplyr %>% left_join filter rename
#' @importFrom tidyr crossing
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
FindStableMarkers <- function(
    id1,
    id2,
    selection_frequencies,
    PFER_max,
    assumption="none",
    B=NULL,
    theta=NULL
){
  if(!is.na(id2) && id1 == id2) stop("id1 must be different from id2!")

  selection_frequencies_filtered <- selection_frequencies %>%
    {
      if(is.na(id2)){
        filter(., id1 == !!id1, is.na(id2))
      } else {
        filter(., (id1 == !!id1 & id2 == !!id2) | (id1 == !!id2 & id2 == !!id1))
      }
    }

  if (!getDoParRegistered()) registerDoSEQ()

  selection_frequency_threshold <- foreach(
    q=unique(selection_frequencies_filtered$q),
    .combine="bind_rows"
  ) %dorng% {
    CalculateSelectionFrequencyThreshold(
      p=length(unique(selection_frequencies$marker)),
      q=q,
      PFER_max=PFER_max,
      B=B,
      assumption=assumption,
      theta=theta
    )
  }

  # join threshold to selection freqs and filter for stable markers
  df_stable_markers <- selection_frequency_threshold %>%
    left_join(
      selection_frequencies_filtered,
      by="q",
      relationship="many-to-many"
    ) %>%
    filter(selection_freq >= selection_freq_thresh)

  # ensure all settings are present (NA if no marker has been selected)
  if(all(is.na(selection_frequencies$p_val_thresh))){
    # fixed q per subsample was selected
    result <- left_join(
      crossing(
        id1=as(id1, class(df_stable_markers$id1)),
        id2=as(id2, class(df_stable_markers$id2)),
        PFER_max=PFER_max,
        assumption=assumption,
        q=unique(selection_frequencies_filtered$q)
      ) %>%
        { # needed if target class of id is factor (otherwise all levels are present)
          if(is.na(id2)){
            filter(., id1 == !!id1, is.na(id2))
          } else {
            filter(., (id1 == !!id1 & id2 == !!id2) | (id1 == !!id2 & id2 == !!id1))
          }
        },
      df_stable_markers,
      by=c("id1", "id2", "PFER_max", "assumption", "q"),
      relationship="one-to-many"
    )
  } else {
    # q is mean of number of selected variables in subsamples based on p_val_thresh
    result <- left_join(
      crossing(
        id1=as(id1, class(df_stable_markers$id1)),
        id2=as(id2, class(df_stable_markers$id2)),
        PFER_max=PFER_max,
        assumption=assumption,
        p_val_thresh=unique(selection_frequencies_filtered$p_val_thresh)
      ) %>%
        { # needed if target class of id is factor (otherwise all levels are present)
          if(is.na(id2)){
            filter(., id1 == !!id1, is.na(id2))
          } else {
            filter(., (id1 == !!id1 & id2 == !!id2) | (id1 == !!id2 & id2 == !!id1))
          }
        },
      df_stable_markers,
      by=c("id1", "id2", "PFER_max", "assumption", "p_val_thresh"),
      relationship="one-to-many"
    )
  }

  # for global markers delete column id2 and rename id1 to id
  if(is.na(id2)){
    result <- result %>%
      select(-id2) %>%
      rename(id=id1)
  }

  return(result)
}



#' Find Global Stable Markers (all clusters)
#'
#' Given selection frequencies and a corresponding threshold, this function finds
#' globally (1 vs rest) stable markers for all clusters. Parallelization can be
#' achieved by setting a futures plan, e.g., \code{plan(multicore)}.
#'
#' @param selection_frequencies A data frame with a columns for the \eqn{p}-value
#'        \eqn{\lambda}, cluster ID1, cluster ID2 (\code{NA} means global marker
#'        for ID1), and columns corresponding to the gene names in the data set.
#'        Values are the selection frequencies.
#' @param PFER_max A numeric value for the maximum tolerable per-family error
#'        rate (PFER).
#' @param assumption A character vector specifying which assumption is made
#'        when calculating the threshold. Options are \code{"none"} for the
#'        Meinshausen & Bühlmann bound and \code{"unimodal"} & \code{"r_concavity"}
#'        for the two different bounds derived by Shah & Samworth. Defaults to
#'        \code{"none"}.
#' @param B The number of splits to generate pairs of complementary subsamples.
#'        Defaults to \code{NULL}, which is only permissible when \code{assumption == "none"}.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability. Defaults to
#'        \code{NULL}, which is only permissible when \code{assumption == "none"}.
#'        For other assumptions the ratio \eqn{q/p} will be used.
#'
#' @return A data frame with columns for cluster ID, \eqn{p}-value threshold,
#'         selection frequency threshold, marker, and observed selection frequency.
#'         \code{NA} entries mean no marker was detected for this combination of
#'         cluster, \eqn{p}-value threshold and selection frequency threshold.
#'
#' @seealso \code{\link{FindStableMarkers}}, \code{\link{CalculateSelectionFrequencyThreshold}}
#'
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
FindAllGlobalStableMarkers <- function(
    selection_frequencies,
    PFER_max,
    assumption="none",
    B=NULL,
    theta=NULL
){
  unique_ids <- unique(selection_frequencies$id1)

  if (!getDoParRegistered()) registerDoSEQ()
  result <- foreach(
    id=unique_ids,
    .combine="rbind",
    .packages=c("dplyr", "reshape2", "tidyr"),
    .export="FindStableMarkers"
  ) %dorng% {
      FindStableMarkers(
        id1=id,
        id2=NA,
        selection_frequencies=selection_frequencies,
        PFER_max=PFER_max,
        B=B,
        assumption=assumption,
        theta=theta
      )
    }

  return(result)
}


#' Find All Local Stable Markers (all pairs of clusters)
#'
#' Given selection frequencies and a corresponding threshold, this function finds
#' locally (1 vs 1) stable markers for all pairs of clusters. Parallelization can
#' be achieved by setting a futures plan, e.g., \code{plan(multicore)}.
#'
#' @param selection_frequencies A data frame with a columns for the \eqn{p}-value
#'        \eqn{\lambda}, cluster ID1, cluster ID2 (\code{NA} means global marker
#'        for ID1), and columns corresponding to the gene names in the data set.
#'        Values are the selection frequencies.
#' @param PFER_max A numeric value for the maximum tolerable per-family error
#'        rate (PFER).
#' @param assumption A character vector specifying which assumption is made
#'        when calculating the threshold. Options are \code{"none"} for the
#'        Meinshausen & Bühlmann bound and \code{"unimodal"} & \code{"r_concavity"}
#'        for the two different bounds derived by Shah & Samworth. Defaults to
#'        \code{"none"}.
#' @param B The number of splits to generate pairs of complementary subsamples.
#'        Defaults to \code{NULL}, which is only permissible when \code{assumption == "none"}.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability. Defaults to
#'        \code{NULL}, which is only permissible when \code{assumption == "none"}.
#'        For other assumptions the ratio \eqn{q/p} will be used.
#'
#' @return A data frame with columns for cluster IDs, \eqn{p}-value threshold,
#'         selection frequency threshold, marker, and observed selection frequency.
#'         \code{NA} entries mean no marker was detected for this combination of
#'         clusters, \eqn{p}-value threshold and selection frequency threshold.
#'
#' @seealso \code{\link{FindStableMarkers}}, \code{\link{CalculateSelectionFrequencyThreshold}}
#'
#' @importFrom utils combn
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
FindAllLocalStableMarkers <- function(
    selection_frequencies,
    PFER_max,
    assumption="none",
    B=B,
    theta=NULL
){
  unique_ids = unique(selection_frequencies$id1)
  n_ids <- length(unique_ids)
  pair_idx <- combn(n_ids, 2)
  settings <- data.frame(
    id1=unique_ids[pair_idx[1, ]],
    id2=unique_ids[pair_idx[2, ]]
  )

  if (!getDoParRegistered()) registerDoSEQ()
  result <- foreach(
    i=1:nrow(settings),
    .combine="rbind",
    .packages=c("dplyr", "reshape2", "tidyr"),
    .export="FindStableMarkers"
  ) %dorng% {
    FindStableMarkers(
      id1=settings$id1[i],
      id2=settings$id2[i],
      selection_frequencies=selection_frequencies,
      PFER_max=PFER_max,
      assumption=assumption,
      B=B,
      theta=theta
    )
  }

  return(result)
}


