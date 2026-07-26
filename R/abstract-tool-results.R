#' Base class for parsed external-tool results
#'
#' Every result exposes `$annotate(idata)`. Results that cannot annotate
#' `ImmunData` override the default with a result-specific explanation and next
#' step; annotatable results override it with their annotation implementation.
#'
#' @keywords internal
ToolResults <- R6::R6Class(
  classname = "ToolResults",
  cloneable = FALSE,
  public = list(
    initialize = function(run) {
      checkmate::assert_r6(run, classes = "ToolRun")
      private$.run <- run
    },
    annotate = function(idata, ...) {
      checkmate::assert_r6(idata, classes = "ImmunData")
      private$annotate_impl(idata, ...)
    }
  ),
  active = list(
    run = function(value) {
      if (!missing(value)) {
        rlang::abort("`run` is read-only.")
      }
      private$.run
    },
    tool = function(value) {
      if (!missing(value)) {
        rlang::abort("`tool` is read-only.")
      }
      private$.run$tool
    }
  ),
  private = list(
    .run = NULL,
    annotate_impl = function(idata, ...) {
      rlang::abort(
        paste(
          sprintf(
            "These `%s` results cannot annotate `ImmunData`.",
            private$.run$tool
          ),
          "The result does not provide an ImmunData annotation.",
          "See the result documentation for available outputs and next steps."
        ),
        class = "immundapter_unannotatable_results"
      )
    }
  )
)
