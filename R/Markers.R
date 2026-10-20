#' Perform Marker Detection Testing in Subsamples
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
#'        will be passed to \code{marker_fn()} as weights.
#' @param marker_fn A function which calculates \eqn{p}-values for a subsample of
#'        the dataset. It is required to have the following correctly named
#'        arguments:
#'        \itemize{
#'          \item \code{data}: The subsample of the dataset.
#'          \item \code{identities}: A vector of cluster labels of the subsample.
#'          \item \code{weights}: The cell stability scores will be handed over
#'                for weighted testing. Can be ignored in the function body if
#'                unweighted testing is desired.
#'          \item \code{cluster_id1}: A single cluster label from in \code{identities}
#'                for which the test is performed.
#'          \item \code{cluster_id2}: A single cluster label from in \code{identities}
#'                for which the test is performed. If \code{NA} is used, a global
#'                (1 vs rest) test for cluster \code{test_id1} should be performed,
#'                otherwise a local (1 vs 1) test for clusters \code{cluster_id1}
#'                and \code{cluster_id2} should be performed.
#'          \item \code{marker_args}: Additional arguments that may be used.
#'        }
#'        It should return a named vector of \eqn{p}-values with the corresponding
#'        gene as name.
#' @param marker_args Additional objects that will be passed to \code{marker_fn}.
#'        Defaults to \code{NULL}.
#' @param marker_global A boolean flag whether to perform testing for global markers
#'        (1 vs rest).
#' @param marker_local A boolean flag whether to perform testing for local markers
#'        (1 vs 1).
#' @param cells_in_columns A boolean flag indicating whether the columns
#'        of \code{data} correspond to observations (i.e. cells) or not. Defaults
#'        to \code{TRUE}.
#' @param packages A character vector of package names that are used by
#'        \code{marker_fn}. Needed for parallelization with \code{doParallel}.
#'        Common choice may be \code{"Seurat"}. Defaults to \code{NULL}.
#' @param export A character vector of (custom) object names that are used by
#'        \code{marker_fn} and are not covered by \code{packages}. Needed for
#'        parallelization with \code{doParallel}. Defaults to \code{NULL}.
#' @param verbose A boolean flag indication whether or not a status update should
#'        be given. Defaults to \code{TRUE}.
#'
#' @return A data frame with \eqn{p}- and log-fold change values. columns are
#'         cluster ID1, cluster ID2 (\code{NA} indicates global marker testing procedure),
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
SubsampleMarkers <- function(
    data,
    reference_identities,
    subsampling_identities,
    marker_alignment,
    cell_stability_score,
    marker_fn,
    marker_args=NULL,
    marker_global=TRUE,
    marker_local=FALSE,
    cells_in_columns=TRUE,
    packages=NULL,
    export=NULL,
    verbose=TRUE
){
  if(!marker_global & !marker_local){ stop("At least of `marker_global` and `marker_local` must be true!") }

  if(!cells_in_columns){
    data <- t(data)
  }

  unique_reference_identities <- sort(unique(reference_identities))


  if (!getDoParRegistered()) registerDoSEQ()

  data_lognorm <- LogNormalizeData(data)

  # test local markers
  if(marker_local){
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

      p_val <- marker_fn(
        data=if(cells_in_columns) subsample else t(subsample),
        identities=subsampling_identities_k,
        weights=cell_stability_k,
        cluster_id1=subsampling_cluster_id1,
        cluster_id2=subsampling_cluster_id2,
        marker_args=marker_args
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
  if(marker_global){
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



      p_val <- marker_fn(
        data=if(cells_in_columns) subsample else t(subsample),
        identities=subsampling_identities_k,
        weights=cell_stability_k,
        cluster_id1=subsampling_cluster_id1,
        cluster_id2=NULL,
        marker_args=marker_args
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



#' Find Stable Markers
#'
#' Given marker selection frequencies and a corresponding threshold, this function finds
#' globally (1 vs rest) or locally (1 vs 1) stable markers.
#'
#' @param id1,id2 Cluster IDs for which the markers are detected. If \code{id2}
#'        is NULL, globally stable markers are found.
#' @param marker_selection_frequencies A data frame with a columns for the \eqn{p}-value
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
    marker_selection_frequencies,
    PFER_max,
    assumption="none",
    B=NULL,
    theta=NULL
){
  if(!is.null(id2) && id1 == id2) stop("id1 must be different from id2!")

  selection_frequencies_filtered <- marker_selection_frequencies %>%
    {
      if(is.null(id2)){
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
      p=length(unique(marker_selection_frequencies$marker)),
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
  if(all(is.na(marker_selection_frequencies$p_val_thresh))){
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
          if(is.null(id2)){
            filter(., id1 == !!id1, is.null(id2))
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
        id2=if(is.null(id2)) NA else as(id2, class(df_stable_markers$id2)),
        PFER_max=PFER_max,
        assumption=assumption,
        p_val_thresh=unique(selection_frequencies_filtered$p_val_thresh)
      ) %>%
        { # needed if target class of id is factor (otherwise all levels are present)
          if(is.null(id2)){
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
  if(is.null(id2)){
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
#' @param marker_selection_frequencies A data frame with a columns for the \eqn{p}-value
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
    marker_selection_frequencies,
    PFER_max,
    assumption="none",
    B=NULL,
    theta=NULL
){
  unique_ids <- unique(marker_selection_frequencies$id1)

  if (!getDoParRegistered()) registerDoSEQ()
  result <- foreach(
    id=unique_ids,
    .combine="rbind",
    .packages=c("dplyr", "tidyr"),
    .export="FindStableMarkers"
  ) %dorng% {
      FindStableMarkers(
        id1=id,
        id2=NULL,
        marker_selection_frequencies=marker_selection_frequencies,
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
#' @param marker_selection_frequencies A data frame with a columns for the \eqn{p}-value
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
    marker_selection_frequencies,
    PFER_max,
    assumption="none",
    B=NULL,
    theta=NULL
){
  unique_ids = unique(marker_selection_frequencies$id1)
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
    .packages=c("dplyr", "tidyr"),
    .export="FindStableMarkers"
  ) %dorng% {
    FindStableMarkers(
      id1=settings$id1[i],
      id2=settings$id2[i],
      marker_selection_frequencies=marker_selection_frequencies,
      PFER_max=PFER_max,
      assumption=assumption,
      B=B,
      theta=theta
    )
  }

  return(result)
}


