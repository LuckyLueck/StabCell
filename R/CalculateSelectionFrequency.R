#' Calculate Selection Frequency
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
CalculateSelectionFrequency <- function(
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


