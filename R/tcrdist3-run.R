#' Completed tcrdist3 execution
#'
#' `Tcrdist3Run` is returned by [Tcrdist3Tool]`$run()`. It contains immutable
#' execution metadata. `$results()` returns [Tcrdist3ClusterResults] or
#' [Tcrdist3DistanceResults], according to the run's configured mode.
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
