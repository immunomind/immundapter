#' Completed RAxML-NG execution
#'
#' `RaxmlNgRun` records one completed family-level analysis: the settings used,
#' an input summary, the output directory, and paths to generated files. It is
#' not a tree. Call `$results()` to obtain [RaxmlNgResults].
#'
#' @keywords internal
RaxmlNgRun <- R6::R6Class(
  classname = "RaxmlNgRun",
  inherit = ToolRun,
  cloneable = FALSE,
  private = list(
    results_impl = function() {
      RaxmlNgResults$new(self)
    }
  )
)
