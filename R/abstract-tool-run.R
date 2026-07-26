#' Immutable record of one external-tool execution
#'
#' A `ToolRun` owns resolved execution metadata and artifact paths. Concrete run
#' classes implement `$results()` so results can be recovered without retaining
#' or recreating the tool object that launched the run.
#'
#' @field status Read-only execution status: `"success"` or `"partial"`.
#' @field tool Read-only name of the tool that created the run.
#' @field path Read-only path to the run directory.
#' @field artifacts Read-only named list of result artifact paths.
#' @field config Read-only plain-list snapshot of the configuration used for
#'   the run.
#' @field input Read-only provenance summary of the input.
#'
#' @keywords internal
ToolRun <- R6::R6Class(
  classname = "ToolRun",
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an immutable record of a tool execution.
    #' @param values A named list with `status`, `tool`, `path`, `artifacts`,
    #'   `config`, and `input` entries.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_choice(values$status, c("success", "partial"))
      checkmate::assert_string(values$tool)
      checkmate::assert_string(values$path)
      checkmate::assert_list(values$artifacts, names = "unique")
      checkmate::assert_list(values$config, names = "unique")
      checkmate::assert_list(values$input, names = "unique")

      private$.status <- values$status
      private$.tool <- values$tool
      private$.path <- values$path
      private$.artifacts <- values$artifacts
      private$.config <- values$config
      private$.input <- values$input
    },
    #' @description
    #' Read the artifacts as a result object. Concrete run classes select and
    #' create the appropriate result class.
    #' @return A [ToolResults] object for this run.
    results = function() {
      private$results_impl()
    }
  ),
  active = list(
    status = function(value) {
      if (!missing(value)) {
        rlang::abort("`status` is read-only.")
      }
      private$.status
    },
    tool = function(value) {
      if (!missing(value)) {
        rlang::abort("`tool` is read-only.")
      }
      private$.tool
    },
    path = function(value) {
      if (!missing(value)) {
        rlang::abort("`path` is read-only.")
      }
      private$.path
    },
    artifacts = function(value) {
      if (!missing(value)) {
        rlang::abort("`artifacts` is read-only.")
      }
      private$.artifacts
    },
    config = function(value) {
      if (!missing(value)) {
        rlang::abort("`config` is read-only.")
      }
      private$.config
    },
    input = function(value) {
      if (!missing(value)) {
        rlang::abort("`input` is read-only.")
      }
      private$.input
    }
  ),
  private = list(
    .status = NULL,
    .tool = NULL,
    .path = NULL,
    .artifacts = NULL,
    .config = NULL,
    .input = NULL,
    results_impl = function() {
      rlang::abort(sprintf(
        "`%s` must implement `$results()` for this run type.",
        class(self)[[1L]]
      ))
    }
  )
)
