#' Completed ImmGLIPH execution
#'
#' `ImmGliphRun` records the settings, input summary, output directory, and
#' generated files for one ImmGLIPH analysis. Call `$results()` to obtain
#' [ImmGliphResults].
#'
#' @keywords internal
ImmGliphRun <- R6::R6Class(
  classname = "ImmGliphRun",
  inherit = ToolRun,
  cloneable = FALSE,
  private = list(
    results_impl = function() {
      ImmGliphResults$new(self)
    }
  )
)
