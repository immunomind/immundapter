#' Completed IgBLAST execution
#'
#' [IgblastTool]`$run()` returns an `IgblastRun`. It records the settings and
#' where the files were saved. Call `$results()` to obtain [IgblastResults].
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
