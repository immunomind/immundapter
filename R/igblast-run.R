#' Completed IgBLAST execution
#'
#' `IgblastRun` is returned by [IgblastTool]`$run()`. It contains immutable
#' execution metadata and artifact paths. `$results()` returns
#' [IgblastResults].
#'
#' @keywords internal
IgblastRun <- R6::R6Class(
  classname = "IgblastRun",
  inherit = ToolRun,
  cloneable = FALSE,
  private = list(
    results_impl = function() {
      IgblastResults$new(self)
    }
  )
)
