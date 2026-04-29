#' Run StabCell Step
#'
#' Runs a single step from StabCell. If the step was already run before, i.e.,
#' the results are already available in the provided \code{checkpoint}, the
#' corresponding results are returned.
#'
#' @param checkpoint List with step names as keys and corresponding results as
#'        values.
#' @param step_index Numeric value of current step.
#' @param steps_to_run Numeric vector of all steps to be run.
#' @param step_name Character value of the current step name.
#' @param step_fn Callable function which returns the results of a current step.
#'
#' @return Result from single StabCell step.
#'
#' @keywords internal
run_step <- function(checkpoint, step_index, steps_to_run, step_name, step_fn) {
  # Check if the current step is in steps_to_run
  if (
    !(step_index %in% steps_to_run)
    || (step_name %in% names(checkpoint))
  ) {
    # returns NULL if "step_index" not in "steps_to_run" and "step_name" not in "checkpoint"
    return(checkpoint[[step_name]])
  } else {
      return(step_fn())
  }
}

#' Validates a StabCell Checkpoint
#'
#' Iterates over all steps in \code{steps_to_run} and checks if the required input
#' values for the variables are available in the provided \code{checkpoint} or if
#' they are created by the other steps. If required inputs are missing an error
#' will be raised.
#'
#' @param steps_to_run Numeric vector of all steps to be run.
#' @param checkpoint List with step names as keys and corresponding results as
#'        values.
#'
#' @keywords internal
validate_checkpoint <- function(steps_to_run, checkpoint){
  required_variables <- character()
  created_variables <- character()
  for(step_index in steps_to_run){
    switch(step_index,
           { # 1
             created_variables <- c(created_variables, "subsampling_identities")
           },
           { # 2
             created_variables <- c(created_variables, "clustering_alignment")
             required_variables <- c(required_variables, "subsampling_identities")
           },
           { # 3
             created_variables <- c(created_variables, "assignment_probabilities")
             required_variables <- c(required_variables, "subsampling_identities", "clustering_alignment")
           },
           { # 4
             created_variables <- c(created_variables, "clustering_mode")
             required_variables <- c(required_variables, "assignment_probabilities")
           },
           { # 5
             created_variables <- c(created_variables, "cell_stability_score")
             required_variables <- c(required_variables, "assignment_probabilities")
           },
           { # 6
             created_variables <- c(created_variables, "marker_alignment")
             required_variables <- c(required_variables, "clustering_mode", "subsampling_identities")
           },
           { # 7
             created_variables <- c(created_variables, "subsampling_p_values")
             required_variables <- c(required_variables, "clustering_mode", "subsampling_identities",
                                     "marker_alignment", "cell_stability_score")
           },
           { # 8
             created_variables <- c(created_variables, "selection_frequencies")
             required_variables <- c(required_variables, "subsampling_p_values")
           },
           { # 9
             created_variables <- c(created_variables, "stable_markers_global")
             required_variables <- c(required_variables, "selection_frequencies")
           },
           { # 10
             created_variables <- c(created_variables, "stable_markers_local")
             required_variables <- c(required_variables, "selection_frequencies")
           }
    )
  }

  required_variables <- setdiff(unique(required_variables), created_variables)
  missing_variables <- setdiff(required_variables, names(checkpoint))

  if(length(missing_variables) > 0){
    stop(paste0(
      "The following intermediate results are missing from 'checkpoint' to run the specified steps: ",
      paste(missing_variables, collapse=", ")
    ))
  }
  return(invisible())
}

#' StabCell
#'
#' This function performs an integrative workflow of clustering and detection of
#' differentially expressed marker genes. Performed steps are
#' \enumerate{
#'    \item \code{\link{SubsampleClustering}}
#'    \item \code{\link{ClusterStabilityAlignment}}
#'    \item \code{\link{CalculateAssignmentProbabilities}}
#'    \item \code{\link{GetSubsampleClusteringMode}}
#'    \item \code{\link{CalculateCellStabilityScore}}
#'    \item \code{\link{MarkerStabilityAlignment}}
#'    \item \code{\link{SubsampleTesting}}
#'    \item \code{\link{CalculateSelectionFrequency}}
#'    \item \code{\link{FindAllGlobalStableMarkers}}
#'    \item \code{\link{FindAllLocalStableMarkers}}
#' }
#' See References for more details. Parallelization can be achieved by setting a
#' futures plan, e.g., \code{plan(multicore)} and calling \code{doFuture::registerDoFuture()}
#' (requires installation of package \code{doFuture}).
#'
#' @param data A matrix- or data frame-like object containing the counts of
#'        gene-cell-pairs.
#' @param reference_identities A vector of initial cluster labels. Usually a
#'        clustering of observations independent of \emph{StabCell} is used here.
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
#'        each observation in the subsample. Factor labels will be disregarded for
#'        their numeric level indicator.
#' @param test_fn A function which calculates \eqn{p}-values for a subsample of
#'        the dataset. It is required to have the following correctly named
#'        arguments:
#'        \itemize{
#'          \item \code{data}: The subsample of the dataset.
#'          \item \code{identities}: A vector of cluster labels of the subsample.
#'          \item \code{weights}: The observation stability scores will be
#'                handed over for weighted testing. Can be ignored in the function
#'                body if unweighted testing is desired.
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
#' @param cluster_args Additional objects that will be passed to \code{cluster_fn}.
#'        Defaults to \code{NULL}.
#' @param test_args Additional objects that will be passed to \code{test_fn}.
#'        Defaults to \code{NULL}.
#' @param test_global A boolean flag whether to perform testing for global markers
#'        (1 vs rest). Defaults to \code{TRUE}.
#' @param test_local A boolean flag whether to perform testing for local markers
#'        (1 vs 1). Defaults to \code{FALSE}.
#' @param PFER_max The maximum number of expected variables with low selection
#'        probability being selected.
#' @param p_value_threshold A numeric value indicating which genes are treated as
#'        markers in each subsample's testing procedure. Defaults to \code{0.05}.
#' @param n_subsample_markers How many variables should be selected in each
#'        subsample. Should be large enough to select all informative marker genes.
#' @param subsampling_iterations An integer specifying how often the subsampling
#'        procedure is performed. Defaults to 50.
#' @param subsampling_fraction A decimal number between 0 and 1 indicating which
#'        fraction of the data will be analyzed in the subsamples. Defaults to 0.5.
#' @param cell_stability_metric A metric how to calculate the stability scores.
#'        Options are \code{"gini"}, \code{"entropy"}, and \code{"L2score"}.
#'        Defaults to \code{"gini"}.
#' @param complementary_pairs A boolean flag indicating whether or not to use
#'        complementary pairs subsampling. If \code{TRUE} (default), twice as
#'        many subsamples as specified by \code{subsampling_iterations} will be
#'        analyzed. Note, that if \code{subsampling_fraction} is not equal to
#'        0.5, the subsamples will be different in size: One will contain
#'        \code{subsampling_fraction} % and the other \code{1 - subsampling_fraction}
#'        of the data.
#' @param assumption A character vector specifying which assumption(s) are made
#'        when calculating the threshold. Options are \code{"none"} for the
#'        Meinshausen & Bühlmann bound and \code{"unimodal"} & \code{"r_concavity"}
#'        for the two different bounds derived by Shah & Samworth. Defaults to
#'        \code{"none"}.
#' @param ranking How to rank the genes. The top \code{n_subsample_markers} will
#'        be selected for further analysis. Options are \code{"p_val"} (increasing
#'        \eqn{p}-value) and \code{"lfc"} (dereasing log-fold change). The other
#'        one is used in case of ties. Defaults to \code{"p_val"}.
#' @param theta The probability threshold to divide the variables into the set of
#'        variables with low and with high selection probability. Defaults to
#'        \code{NULL}, which is only permissible is \code{assumption == "none"}.
#'        For other assumptions the ratio \eqn{q/p} will be used.
#' @param cells_in_columns A boolean flag indicating whether the columns
#'        of \code{data} correspond to observations (i.e. cells) or not. Defaults
#'        to \code{TRUE}.
#' @param packages A character vector of package names that are used by
#'        \code{cluster_fn} and \code{test_fn}. Needed for parallelization worker
#'        setup. Common choice may be \code{"Seurat"}. Defaults to \code{NULL}.
#' @param export A character vector of (custom) object names that are used by
#'        \code{cluster_fn} & \code{test_fn} and are not covered by
#'        \code{packages}. Needed for parallelization worker setup. Defaults to
#'        \code{NULL}.
#' @param steps_to_run A vector of integers indicating which of the analysis steps
#'        should be performed (see Description). If step 1 is not included, or
#'        any other intermediate step is skipped, the corresponding step results
#'        should be included in the \code{checkpoint} argument. Defaults to
#'        \code{NULL}, i.e., every step is performed.
#' @param checkpoint A named list of intermediary results (see "Value" section).
#'        Existing values of names will be used for further calculations, while
#'        missing names will be added. Defaults to \code{NULL}.
#' @param checkpoint_path A path to a \code{RData}-file, in which intermediary
#'        results will be saved. Defaults to \code{NULL}.
#' @param verbose A boolean flag whether to give occasional updates which step
#'        is currently performed. Defaults to \code{TRUE}.
#'
#' @return A named list with results of each step that is performed
#'         \itemize{
#'            \item \code{subsampling_identities}: A matrix of cluster labels for
#'                  each subsample.
#'            \item \code{clustering_alignment}: A data frame which aligns a
#'                  cluster from \code{reference_identities} to every cluster in
#'                  every subsample clustering.
#'            \item \code{assignment_probabilities}: A data frame with assignment
#'                  probabilities for every observation to every cluster in
#'                  \code{reference_identities}.
#'            \item \code{clustering_mode}: A vector of cluster labels
#'                  corresponding to the highest \code{assignment_probabilities}.
#'            \item \code{cell_stability_score}: A numeric vector with
#'                  stability scores of every cell. Higher scores indicate more
#'                  consistent assignment to the same reference cluster.
#'            \item \code{marker_alignment}: A data frame which aligns a cluster
#'                  from every subsample to every cluster in \code{clustering_mode}.
#'            \item \code{subsampling_p_values}: A data frame with \eqn{p}-values
#'                  for every cluster (pair) in \code{clustering_mode} and every
#'                  gene in every subsample.
#'            \item \code{selection_frequencies}: A data frame with selection
#'                  frequencies for every gene and every cluster (pair) depending
#'                  on \code{p_value_threshold}.
#'            \item \code{stable_markers_global}: A data frame with globally
#'                  stable markers for each cluster in \code{clustering_mode}
#'                  depending on \code{p_value_threshold} and
#'                  \code{selection_frequency_threshold}.
#'            \item \code{stable_markers_local}: A data frame with locally stable
#'                  markers for each pair of clusters in \code{clustering_mode}
#'                  depending on \code{p_value_threshold} and
#'                  \code{selection_frequency_threshold}.
#'         }
#'
#' @export
StabCell <- function(
    data,
    reference_identities,
    cluster_fn,
    test_fn,
    PFER_max,
    p_value_threshold=0.05,
    n_subsample_markers=NULL,
    cluster_args=NULL,
    test_args=NULL,
    test_global=TRUE,
    test_local=FALSE,
    subsampling_iterations=50,
    subsampling_fraction=0.5,
    complementary_pairs=TRUE,
    cell_stability_metric=c("gini", "entropy", "L2score"),
    assumption="none",
    ranking="p_val",
    theta=NULL,
    cells_in_columns=TRUE,
    packages=NULL,
    export=NULL,
    steps_to_run=NULL,
    checkpoint=NULL,
    checkpoint_path=NULL,
    verbose=TRUE
){
  if(!complementary_pairs && length(intersect(assumption, c("unimodal", "r_concavity"))) > 0){
    stop("Unimodal and r-concavity assumptions require complementary pairs!")
  }

  # only perform subset of steps?
  if(is.null(steps_to_run)){
    steps_to_run <- 1:10
  }

  if(is.null(checkpoint)){
    # start from empty checkpoint
    checkpoint <- list()
  } else {
    # check if checkpoint is valid for specified steps
    validate_checkpoint(checkpoint=checkpoint, steps_to_run=steps_to_run)
  }


  # Clustering
  subsampling_identities <- run_step(
    checkpoint=checkpoint,
    step_index=1,
    steps_to_run=steps_to_run,
    step_name="subsampling_identities",
    step_fn=function(){
      if(verbose) message("subsampling & clustering...\n")
      subsampling_identities <- SubsampleClustering(
        data=data,
        cluster_fn=cluster_fn,
        cluster_args=cluster_args,
        subsampling_fraction=subsampling_fraction,
        subsampling_iterations=subsampling_iterations,
        complementary_pairs=complementary_pairs,
        observations_in_columns=cells_in_columns,
        packages=packages,
        export=export
      )
      checkpoint[["subsampling_identities"]] <<- subsampling_identities

      if(!is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(subsampling_identities)
    }
  )

  # Clustering Stability
  clustering_alignment <- run_step(
    checkpoint=checkpoint,
    step_index=2,
    steps_to_run=steps_to_run,
    step_name="clustering_alignment",
    step_fn=function(){
      if(verbose) message("clustering stability...\n")
      clustering_alignment <- ClusterStabilityAlignment(
        reference_identities=reference_identities,
        subsampling_identities=subsampling_identities
      )
      checkpoint[["clustering_alignment"]] <<- clustering_alignment
      return(clustering_alignment)
    }
  )

  assignment_probabilities <- run_step(
    checkpoint=checkpoint,
    step_index=3,
    steps_to_run=steps_to_run,
    step_name="assignment_probabilities",
    step_fn=function(){
      assignment_probabilities <- CalculateAssignmentProbabilities(
        subsampling_identities=subsampling_identities,
        clustering_alignment=clustering_alignment
      )
      checkpoint[["assignment_probabilities"]] <<- assignment_probabilities
      return(assignment_probabilities)
    }
  )

  clustering_mode <- run_step(
    checkpoint=checkpoint,
    step_index=4,
    steps_to_run=steps_to_run,
    step_name="clustering_mode",
    step_fn=function(){
      clustering_mode <- GetSubsampleClusteringMode(assignment_probabilities)
      checkpoint[["clustering_mode"]] <<- clustering_mode
      if(!is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(clustering_mode)
    }
  )

  cell_stability_score <- run_step(
    checkpoint=checkpoint,
    step_index=5,
    steps_to_run=steps_to_run,
    step_name="cell_stability_score",
    step_fn=function(){
      cell_stability_score <- CalculateCellStabilityScore(
        assignment_probabilities=assignment_probabilities,
        metric=cell_stability_metric,
        normalize=TRUE
      )
      checkpoint[["cell_stability_score"]] <<- cell_stability_score

      if(!is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(cell_stability_score)
    }
  )


  # Testing
  marker_alignment <- run_step(
    checkpoint=checkpoint,
    step_index=6,
    steps_to_run=steps_to_run,
    step_name="marker_alignment",
    step_fn=function(){
      if(verbose) message("testing...\n")
      marker_alignment <- MarkerStabilityAlignment(
        reference_identities=clustering_mode,
        subsampling_identities=subsampling_identities
      )
      checkpoint[["marker_alignment"]] <<- marker_alignment
      return(marker_alignment)
    }
  )

  subsampling_p_values <- run_step(
    checkpoint=checkpoint,
    step_index=7,
    steps_to_run=steps_to_run,
    step_name="subsampling_p_values",
    step_fn=function(){
      subsampling_p_values <- SubsampleTesting(
        data=data,
        reference_identities=clustering_mode,
        subsampling_identities=subsampling_identities,
        marker_alignment=marker_alignment,
        cell_stability_score=cell_stability_score,
        test_fn=test_fn,
        test_args=test_args,
        test_global=test_global,
        test_local=test_local,
        observations_in_columns=cells_in_columns,
        packages=packages,
        export=export,
        verbose=verbose
      )
      checkpoint[["subsampling_p_values"]] <<- subsampling_p_values

      if(!is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(subsampling_p_values)
    }
  )

  # selection frequencies
  selection_frequencies <- run_step(
    checkpoint=checkpoint,
    step_index=8,
    steps_to_run=steps_to_run,
    step_name="selection_frequencies",
    step_fn=function(){
      if(verbose) message("calculate selection frequency...\n")
      selection_frequencies <- CalculateSelectionFrequency(
        subsampling_p_values=subsampling_p_values,
        p_value_threshold=p_value_threshold,
        n_subsample_markers=n_subsample_markers,
        ranking=ranking
      )
      checkpoint[["selection_frequencies"]] <<- selection_frequencies

      if(!is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(selection_frequencies)
    }
  )

  # stable markers
  stable_markers_global <- run_step(
    checkpoint=checkpoint,
    step_index=9,
    steps_to_run=steps_to_run,
    step_name="stable_markers_global",
    step_fn=function(){
      if(test_global){
        if(verbose) message("find global markers...\n")
        stable_markers_global <- FindAllGlobalStableMarkers(
          selection_frequencies=selection_frequencies,
          PFER_max=PFER_max,
          assumption=assumption,
          B=subsampling_iterations,
          theta=theta
        )
      } else {
        stable_markers_global <- NULL
      }
      checkpoint[["stable_markers_global"]] <<- stable_markers_global

      if(test_global && !is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(stable_markers_global)
    }
  )

  stable_markers_local <- run_step(
    checkpoint=checkpoint,
    step_index=10,
    steps_to_run=steps_to_run,
    step_name="stable_markers_local",
    step_fn=function(){
      if(test_local){
        if(verbose) message("find local markers...\n")
        stable_markers_local <- FindAllLocalStableMarkers(
          selection_frequencies=selection_frequencies,
          PFER_max=PFER_max,
          assumption=assumption,
          B=subsampling_iterations,
          theta=theta
        )
      } else {
        stable_markers_local <- NULL
      }

      checkpoint[["stable_markers_local"]] <<- stable_markers_local

      if(test_local && !is.null(checkpoint_path)) save(checkpoint, file=checkpoint_path)
      return(stable_markers_local)
    }
  )

  if(verbose) message("done!\n")

  return(checkpoint)
}

