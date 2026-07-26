#' Base class for external-tool adapters
#'
#' `AbstractTool` handles the parts shared by every adapter: configuration,
#' run-directory creation, and immutable `ToolRun` construction. Concrete tools
#' validate and normalize their supported input types, execute the external
#' tool, and select the concrete run class that parses their results.
#'
#' @keywords internal
tool_input_summary <- function(input) {
  if (inherits(input, "ImmunData")) {
    provenance <- input$provenance
    return(list(
      type = "immundata",
      snapshot_id = provenance$snapshot_id,
      snapshot_path = provenance$current_path
    ))
  }

  if (is.character(input)) {
    paths <- as.character(fs::path_abs(input))
    names(paths) <- names(input)
    return(list(
      type = "files",
      paths = paths
    ))
  }

  list(
    type = "object",
    class = class(input)
  )
}

resolve_tool_output_dir <- function(input, output_dir) {
  if (!is.null(output_dir)) {
    return(fs::path_abs(output_dir))
  }

  if (!inherits(input, "ImmunData")) {
    return(fs::path_home())
  }

  artifacts_path <- input$provenance$artifacts_path
  if (is.null(artifacts_path)) {
    rlang::abort(paste(
      "`output_dir` is required when `input` is an `ImmunData` object",
      "without a persisted snapshot artifacts path."
    ))
  }

  fs::path_abs(artifacts_path)
}

create_tool_run_dir <- function(output_dir, tool) {
  fs::dir_create(output_dir)

  stem <- paste0(tool, "-", format(Sys.time(), "%Y%m%d-%H%M%S"))
  run_dir <- fs::path(output_dir, stem)
  suffix <- 0L

  while (fs::dir_exists(run_dir) || fs::file_exists(run_dir)) {
    suffix <- suffix + 1L
    run_dir <- fs::path(
      output_dir,
      paste0(stem, "-", sprintf("%03d", suffix))
    )
  }

  fs::dir_create(run_dir)
  run_dir
}

AbstractTool <- R6::R6Class(
  classname = "AbstractTool",
  cloneable = FALSE,
  public = list(
    initialize = function(name, config, run_class) {
      checkmate::assert_string(name)
      checkmate::assert_r6(config, classes = "AbstractConfig")
      is_tool_run_class <- inherits(run_class, "R6ClassGenerator") &&
        !identical(run_class, ToolRun)
      current_class <- run_class
      while (
        is_tool_run_class &&
          !identical(current_class, ToolRun) &&
          !is.null(current_class$get_inherit())
      ) {
        current_class <- current_class$get_inherit()
      }
      is_tool_run_class <- is_tool_run_class &&
        identical(current_class, ToolRun)
      if (!is_tool_run_class) {
        rlang::abort(
          "`run_class` must be a concrete `ToolRun` R6 class generator."
        )
      }

      private$.name <- name
      private$.config <- config
      private$.config_class <- class(config)[[1L]]
      private$.run_class <- run_class
    },
    help = function() {
      private$help_impl()
    },
    check = function() {
      private$check_impl()
    },
    citation = function() {
      private$citation_impl()
    },
    run = function(input, config = NULL, output_dir = NULL) {
      checkmate::assert_string(output_dir, null.ok = TRUE)

      if (is.null(config)) {
        config <- private$.config
      } else {
        checkmate::assert_r6(config, classes = private$.config_class)
      }

      input <- private$validate_input_impl(input, config)
      output_dir <- resolve_tool_output_dir(input, output_dir)
      run_dir <- create_tool_run_dir(output_dir, private$.name)

      artifacts <- private$run_impl(input, config, run_dir)

      private$.run_class$new(list(
        status = "success",
        tool = private$.name,
        path = run_dir,
        artifacts = artifacts,
        config = config$as_list(),
        input = tool_input_summary(input)
      ))
    }
  ),
  active = list(
    name = function(value) {
      if (!missing(value)) {
        rlang::abort("`name` is read-only.")
      }
      private$.name
    },
    config = function(value) {
      if (!missing(value)) {
        rlang::abort("`config` is read-only.")
      }
      private$.config
    }
  ),
  private = list(
    .name = NULL,
    .config = NULL,
    .config_class = NULL,
    .run_class = NULL,
    help_impl = function() {
      rlang::abort("The tool must implement `help()`.")
    },
    check_impl = function() {
      rlang::abort("The tool must implement `check()`.")
    },
    citation_impl = function() {
      rlang::abort("The tool must implement `citation()`.")
    },
    validate_input_impl = function(input, config) {
      rlang::abort("The tool must implement input validation for `run()`.")
    },
    run_impl = function(input, config, run_dir) {
      rlang::abort("The tool must implement `run()`.")
    }
  )
)
