#' Base class for parsed external-tool results
#'
#' Every result exposes `$annotate(idata)`. Results that cannot annotate
#' `ImmunData` override the default with a result-specific explanation and next
#' step; annotatable results override it with their annotation implementation.
#'
#' @field run The read-only [ToolRun] that produced these results.
#' @field tool The read-only name of the tool that produced these results.
#'
#' @keywords internal
ToolResults <- R6::R6Class(
  classname = "ToolResults",
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create a result object for a completed tool run.
    #' @param run The [ToolRun] that produced the result artifacts.
    initialize = function(run) {
      checkmate::assert_r6(run, classes = "ToolRun")
      private$.run <- run
    },
    #' @description
    #' Add result annotations to immune-repertoire data. Concrete result
    #' classes implement the annotation operation. The base method reports an
    #' informative error when a result cannot provide annotations.
    #' @param idata An [immundata::ImmunData] object to annotate.
    #' @param ... Additional arguments passed to the concrete annotation
    #'   method.
    #' @return The value returned by the concrete annotation method, normally
    #'   an updated [immundata::ImmunData] object.
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
