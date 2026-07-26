#' Base class for external-tool adapters
#'
#' `AbstractTool` handles the parts shared by every adapter: configuration,
#' run-directory creation, the `ToolRun` return value, and result conversion.
#' Concrete tools validate and normalize their supported input types and
#' implement the remaining private hooks.
#'
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
    initialize = function(name, config) {
      checkmate::assert_string(name)
      checkmate::assert_r6(config, classes = "AbstractConfig")

      private$.name <- name
      private$.config <- config
      private$.config_class <- class(config)[[1L]]
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

      structure(
        list(
          status = "success",
          tool = private$.name,
          path = run_dir,
          artifacts = artifacts,
          config = config$as_list(),
          input = tool_input_summary(input)
        ),
        class = c("ToolRun", "list")
      )
    },
    results = function(run, output = c("duckplyr_df", "df")) {
      output <- match.arg(output)
      checkmate::assert_class(run, "ToolRun")

      if (!identical(run$tool, private$.name)) {
        rlang::abort(sprintf(
          "This is a '%s' run, not a '%s' run.",
          run$tool,
          private$.name
        ))
      }

      result <- private$results_impl(run)
      if (identical(output, "df")) {
        return(as.data.frame(dplyr::collect(result)))
      }
      result
    },
    annotate = function(idata, results, ...) {
      checkmate::assert_r6(idata, classes = "ImmunData")
      checkmate::assert_data_frame(results)
      private$annotate_impl(idata, results, ...)
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
    },
    results_impl = function(run) {
      rlang::abort("The tool must implement `results()`.")
    },
    annotate_impl = function(idata, results, ...) {
      rlang::abort("The tool must implement `annotate()`.")
    }
  )
)
