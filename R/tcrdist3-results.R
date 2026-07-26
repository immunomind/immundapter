#' Parsed tcrdist3 results
#'
#' `Tcrdist3Results` contains outputs shared by every tcrdist3 mode: a lazy
#' duckplyr receptor mapping and the path to the NumPy distance matrix.
#'
#' @field mapping Read-only lazy duckplyr table that maps each original input
#'   receptor to the receptor used for distance calculation.
#' @field distances Read-only path to the NumPy `.npy` distance matrix.
#'
#' @keywords internal
Tcrdist3Results <- R6::R6Class(
  classname = "Tcrdist3Results",
  inherit = ToolResults,
  cloneable = FALSE,
  active = list(
    mapping = function(value) {
      if (!missing(value)) {
        rlang::abort("`mapping` is read-only.")
      }
      path <- private$.run$artifacts$mapping
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    },
    distances = function(value) {
      if (!missing(value)) {
        rlang::abort("`distances` is read-only.")
      }
      path <- private$.run$artifacts$distances
      checkmate::assert_file_exists(path)
      path
    }
  )
)


#' Parsed tcrdist3 distance-only results
#'
#' Distance-only results contain a receptor mapping and distance matrix but no
#' receptor-level annotations. `$annotate()` explains how to produce
#' annotatable cluster results.
#'
#' @keywords internal
Tcrdist3DistanceResults <- R6::R6Class(
  classname = "Tcrdist3DistanceResults",
  inherit = Tcrdist3Results,
  cloneable = FALSE,
  private = list(
    annotate_impl = function(idata, ...) {
      rlang::abort(
        paste(
          "These tcrdist3 distance results cannot annotate `ImmunData`",
          "because the run produced a distance matrix and receptor mapping,",
          "but no receptor-level annotations.",
          "Rerun tcrdist3 with `mode = \"cluster\"`, choose an appropriate",
          "`radius`, and call `$annotate()` on those cluster results."
        ),
        class = "immundapter_unannotatable_results"
      )
    }
  )
)


#' Parsed tcrdist3 cluster results
#'
#' Cluster results add lazy duckplyr cluster and edge tables to the common
#' tcrdist3 outputs. `$annotate(idata)` joins receptor-level cluster
#' assignments to an `ImmunData` object as `tcrdist_cluster_id` and
#' `tcrdist_cluster_size`.
#'
#' @field clusters Read-only lazy duckplyr table of receptor cluster
#'   assignments and cluster sizes.
#' @field edges Read-only lazy duckplyr table of receptor pairs connected at
#'   the configured distance radius.
#'
#' @keywords internal
Tcrdist3ClusterResults <- R6::R6Class(
  classname = "Tcrdist3ClusterResults",
  inherit = Tcrdist3Results,
  cloneable = FALSE,
  active = list(
    clusters = function(value) {
      if (!missing(value)) {
        rlang::abort("`clusters` is read-only.")
      }
      path <- private$.run$artifacts$clusters
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    },
    edges = function(value) {
      if (!missing(value)) {
        rlang::abort("`edges` is read-only.")
      }
      path <- private$.run$artifacts$edges
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(path, options = list(delim = "\t"))
    }
  ),
  private = list(
    annotate_impl = function(idata, ...) {
      receptor_col <- immundata::imd_schema("receptor")
      annotations <- self$clusters |>
        dplyr::select(
          dplyr::all_of(c(
            "input_id",
            "tcrdist_cluster_id",
            "tcrdist_cluster_size"
          ))
        ) |>
        dplyr::rename(dplyr::all_of(
          stats::setNames("input_id", receptor_col)
        ))

      immundata::annotate_receptors(
        idata = idata,
        annotations = annotations,
        annot_col = receptor_col,
        conflicts = "replace",
        ...
      )
    }
  )
)
