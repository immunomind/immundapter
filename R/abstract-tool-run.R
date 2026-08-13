#' Record of one external-tool run
#'
#' A `ToolRun` records the settings, input summary, output directory, and files
#' created by one run. Use `$results()` to read the biological results.
#'
#' @field status Read-only execution status: `"success"` or `"partial"`.
#' @field tool Read-only name of the tool that created the run.
#' @field path Read-only path to the run directory.
#' @field artifacts Read-only named list of paths to generated files.
#' @field config Read-only list of settings used for the run.
#' @field input Read-only provenance summary of the input.
#'
#' @keywords internal
ToolRun <- R6::R6Class(
  classname = "ToolRun",
  cloneable = FALSE,
  public = list(
    #' @description
    #' Record a completed tool run.
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
    #' Read the result tables and files.
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
