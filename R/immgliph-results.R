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
      private$artifact("groups", value, type = "parquet")
    },
    connections = function(value) {
      private$artifact("connections", value, type = "parquet")
    },
    cluster_properties = function(value) {
      private$artifact("cluster_properties", value, type = "parquet")
    },
    mapping = function(value) {
      private$artifact("mapping", value, type = "parquet")
    },
    raw_output = function(value) {
      private$artifact("result", value, type = "path", field = "raw_output")
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
