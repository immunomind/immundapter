#' Parsed ImmGLIPH results
#'
#' `ImmGliphResults` provides the receptor groups, connections, group
#' properties, the columns sent to ImmGLIPH, and the original R result file.
#'
#' @field groups Read-only table with one row per grouped receptor.
#' @field connections Read-only table of sequence connections found by
#'   ImmGLIPH.
#' @field cluster_properties Read-only table of group properties and scores.
#' @field mapping Read-only table containing the receptor columns sent to
#'   ImmGLIPH.
#' @field raw_output Read-only path to the original ImmGLIPH `.rds` result.
#'
#' @keywords internal
ImmGliphResults <- R6::R6Class(
  classname = "ImmGliphResults",
  inherit = ToolResults,
  cloneable = FALSE,
  active = list(
    groups = function(value) {
      if (!missing(value)) {
        rlang::abort("`groups` is read-only.")
      }
      path <- private$.run$artifacts$groups
      checkmate::assert_file_exists(path)
      duckplyr::read_parquet_duckdb(path)
    },
    connections = function(value) {
      if (!missing(value)) {
        rlang::abort("`connections` is read-only.")
      }
      path <- private$.run$artifacts$connections
      checkmate::assert_file_exists(path)
      duckplyr::read_parquet_duckdb(path)
    },
    cluster_properties = function(value) {
      if (!missing(value)) {
        rlang::abort("`cluster_properties` is read-only.")
      }
      path <- private$.run$artifacts$cluster_properties
      checkmate::assert_file_exists(path)
      duckplyr::read_parquet_duckdb(path)
    },
    mapping = function(value) {
      if (!missing(value)) {
        rlang::abort("`mapping` is read-only.")
      }
      path <- private$.run$artifacts$mapping
      checkmate::assert_file_exists(path)
      duckplyr::read_parquet_duckdb(path)
    },
    raw_output = function(value) {
      if (!missing(value)) {
        rlang::abort("`raw_output` is read-only.")
      }
      path <- private$.run$artifacts$result
      checkmate::assert_file_exists(path)
      path
    }
  ),
  private = list(
    annotate_impl = function(idata, ...) {
      receptor_col <- immundata::imd_schema("receptor")
      annotations <- self$groups |>
        dplyr::select(dplyr::all_of(c(
          receptor_col,
          "immgliph_group_id"
        )))

      immundata::annotate_receptors(
        idata = idata,
        annotations = annotations,
        annot_col = receptor_col,
        ...
      )
    }
  )
)
