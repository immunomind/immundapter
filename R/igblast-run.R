#' @noRd
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
