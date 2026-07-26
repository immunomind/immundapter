#' Parsed tcrdist3 cluster results
#'
#' `Tcrdist3Results` exposes lazy duckplyr tables for clusters, mappings, and
#' edges, plus the path to the NumPy distance matrix. `$annotate(idata)` joins
#' receptor-level cluster assignments to an `ImmunData` object.
#'
#' @keywords internal
Tcrdist3Results <- R6::R6Class(
  classname = "Tcrdist3Results",
  inherit = ToolResults,
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
    mapping = function(value) {
      if (!missing(value)) {
        rlang::abort("`mapping` is read-only.")
      }
      path <- private$.run$artifacts$mapping
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
    },
    distances = function(value) {
      if (!missing(value)) {
        rlang::abort("`distances` is read-only.")
      }
      path <- private$.run$artifacts$distances
      checkmate::assert_file_exists(path)
      path
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
        dplyr::rename(
          !!receptor_col := dplyr::all_of("input_id")
        )

      immundata::annotate_receptors(
        idata = idata,
        annotations = annotations,
        annot_col = receptor_col,
        ...
      )
    }
  )
)
