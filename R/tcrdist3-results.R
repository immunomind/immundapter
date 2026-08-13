#' Parsed tcrdist3 results
#'
#' `Tcrdist3Results` contains output shared by both analysis modes. `$mapping`
#' shows which identical receptors were calculated once. Use
#' `dplyr::collect()` to read this table into R. `$distances` gives the path to
#' the distance matrix.
#'
#' @field mapping Read-only table showing which original receptors were
#'   identical and calculated once.
#' @field distances Read-only file path to the pairwise distance matrix in
#'   Python NumPy `.npy` format. This field is a path, not an R matrix.
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
#' Distance-only results contain a distance matrix but no receptor groups. Run
#' tcrdist3 again with `mode = "cluster"` before adding groups to `ImmunData`.
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
          "because the run produced distances but no receptor groups.",
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
#' Cluster results add two tables to the common tcrdist3 output.
#' `$clusters` contains the group assignment for each receptor. `$edges`
#' contains pairs of receptors directly connected at the selected distance
#' radius. Use `dplyr::collect()` to bring either table into R memory.
#'
#' `$annotate(idata)` returns a new `ImmunData` dataset with
#' `tcrdist_cluster_id` and `tcrdist_cluster_size` added to receptor
#' columns.
#'
#' @field clusters Read-only table with one cluster number and cluster size per
#'   receptor.
#' @field edges Read-only table of receptor pairs whose distance is within the
#'   selected radius.
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
