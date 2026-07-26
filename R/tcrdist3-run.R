#' Completed tcrdist3 execution
#'
#' `Tcrdist3Run` is returned by [Tcrdist3Tool]`$run()`. It contains immutable
#' execution metadata and creates [Tcrdist3Results] through `$results()`.
#'
#' @keywords internal
Tcrdist3Run <- R6::R6Class(
  classname = "Tcrdist3Run",
  inherit = ToolRun,
  cloneable = FALSE,
  private = list(
    results_impl = function() {
      Tcrdist3Results$new(self)
    }
  )
)
