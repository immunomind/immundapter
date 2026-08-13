#' Completed tcrdist3 execution
#'
#' [Tcrdist3Tool]`$run()` returns a `Tcrdist3Run`. It records the settings and
#' where the files were saved. `$results()` returns [Tcrdist3ClusterResults] or
#' [Tcrdist3DistanceResults], depending on the selected mode.
#'
#' @keywords internal
Tcrdist3Run <- R6::R6Class(
  classname = "Tcrdist3Run",
  inherit = ToolRun,
  cloneable = FALSE,
  private = list(
    results_impl = function() {
      switch(
        self$config$mode,
        cluster = Tcrdist3ClusterResults$new(self),
        distance = Tcrdist3DistanceResults$new(self),
        rlang::abort(sprintf(
          "Unknown tcrdist3 result mode: %s.",
          self$config$mode
        ))
      )
    }
  )
)
