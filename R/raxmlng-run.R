#' Completed RAxML-NG execution
#'
#' `RaxmlNgRun` records the configuration, input provenance, status, and
#' artifacts for one family-level batch. Call `$results()` to obtain
#' [RaxmlNgResults].
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
