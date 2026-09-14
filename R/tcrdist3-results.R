#' @noRd
Tcrdist3Results <- R6::R6Class(
  classname = "Tcrdist3Results",
  inherit = ToolResults,
  cloneable = FALSE,
  active = list(
    mapping = function(value) {
      private$artifact("mapping", value)
    },
    distances = function(value) {
      private$artifact("distances", value, type = "path")
    }
  )
)


#' @noRd
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


#' @noRd
Tcrdist3ClusterResults <- R6::R6Class(
  classname = "Tcrdist3ClusterResults",
  inherit = Tcrdist3Results,
  cloneable = FALSE,
  active = list(
    clusters = function(value) {
      private$artifact("clusters", value)
    },
    edges = function(value) {
      private$artifact("edges", value)
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
