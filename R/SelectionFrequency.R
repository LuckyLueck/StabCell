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
  MAXa <- 1e10
  MINa <- 1e-10
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


#' Calculate Marker Selection Frequency
#'
#' Given \eqn{p}-values for every gene and every cluster in all subsamples,
#' calculate the frequency of \eqn{p}-values being smaller than a specified
#' threshold. One of \code{p_value_threshold} and \code{n_subsample_markers} must
#' be provided. Parallelization can be achieved by setting a futures plan, e.g.,
#' \code{plan(multicore)}.
#'
#' @param subsampling_p_values A \code{data.frame} with a columns for the subsample
#'        ID, cluster ID1, cluster ID2 (\code{NA} means global marker for ID1),
#'        and columns corresponding to the gene names in the data set. Values are
#'        the \eqn{p}-values.
#' @param p_value_threshold In each subsample all genes with a \eqn{p}-value
#'        \code{<= p_value_threshold} are selected as markers. Will overwrite
#'        \code{n_subsample_markers} if not \code{NULL}.
#' @param n_subsample_markers How many genes are selected as markers in every
#'        comparison in every subsample. May be ignored if \code{p_value_threshold}
#'        is not \code{NULL}.
#' @param ranking How to first rank the genes when choosing a fixed number of marker
#'        genes per comparison. Options are \code{"p_val"} (increasing \eqn{p}-value)
#'        and \code{"lfc"} (dereasing log-fold change). The other one is used in
#'        case of ties. Defaults to \code{"p_val"}.
#'
#' @return A \code{data.frame} with a columns for the \eqn{p}-value threshold
#'         \eqn{\lambda}, cluster ID1, cluster ID2 (\code{NA} means global marker
#'         for ID1), marker name, selection frequency, and mean and standard
#'         deviation of \eqn{p}-value and log-fold change.
#'
#' @importFrom stats sd
#' @importFrom dplyr arrange desc summarise n_distinct
#' @importFrom doRNG %dorng%
#' @importFrom foreach foreach getDoParRegistered registerDoSEQ
#'
#' @export
CalculateMarkerSelectionFrequency <- function(
    subsampling_p_values,
    p_value_threshold=NULL,
    n_subsample_markers=NULL,
    ranking=c("p_val", "lfc")
){
  # find selected markers
  if(is.null(p_value_threshold)){
    if(is.null(n_subsample_markers)){
      stop("Either p_value_threshold or n_subsample_markers must be provided!")
    }

    ranking <- match.arg(ranking)

    ### pick top q markers based in every subsample ###
    # -> number of markers is constant between subsamples
    if (!getDoParRegistered()) registerDoSEQ()
    df_selection <- foreach(
      q=n_subsample_markers,
      .combine="bind_rows"
    ) %dorng%{
      temp_df <- subsampling_p_values %>%
        group_by(subsample, id1, id2) %>%
        {
          if(ranking == "p_val"){
            arrange(., p_val, desc(is.finite(lfc)), desc(abs(lfc)))
          } else {
            arrange(., desc(is.finite(lfc)), desc(abs(lfc)), p_val)
          }
        } %>%
        mutate(selected=as.integer(row_number() <= !!q)) %>%   # 1 if in top q, 0 otherwise
        ungroup()

      cbind(p_val_thresh=NA, q=q, temp_df)
    }

  } else {
    if(!is.null(n_subsample_markers)){
      message("Both p_value_threshold and n_subsample_markers were be provided! Proceed only with p-value settings!")
    }
    ### pick top markers based on p value ###
    # -> number of markers may vary between subsamples

    df_selection <- foreach(
      p_thresh=p_value_threshold,
      .combine="bind_rows"
    ) %dorng% {
      subsampling_p_values %>%
        mutate(
          selected=as.integer(p_val <= !!p_thresh),
          p_val_thresh=!!p_thresh
        )
    }

    # calculate expected number of selected markers (i.e. mean)
    df_q <- df_selection %>%
      group_by(p_val_thresh, id1, id2) %>%
      summarise(
        q=sum(selected * 1.0) / n_distinct(subsample),
        .groups="drop"
      ) %>%
      select(p_val_thresh, q, id1, id2)

    df_selection <- left_join(df_q, df_selection, by=c("p_val_thresh", "id1", "id2"))

  }

  # calculate selection frequency per marker
  selection_frequencies <- df_selection %>%
    group_by(p_val_thresh, q, id1, id2, marker) %>%
    summarise(
      selection_freq=sum(selected) / n_distinct(subsample),
      p_val_mean=mean(p_val),
      p_val_sd=sd(p_val),
      lfc_mean=mean(lfc[is.finite(lfc)]),
      lfc_sd=sd(lfc[is.finite(lfc)]),
      .groups="drop"
    )

  return(selection_frequencies)
}

