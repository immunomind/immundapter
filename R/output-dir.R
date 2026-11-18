#' Set the default output directory for file-based runs
#'
#' Save a directory for file-based runs across R sessions. An explicit
#' `output_dir` passed to `$run()` takes precedence. The directory is created
#' when a run first uses it.
#'
#' @param path Directory to save. In an interactive session, omit `path` to
#'   enter and confirm a directory at the console.
#' @return Invisibly, the absolute path that was saved, or `NULL` if cancelled.
#' @export
set_output_dir <- function(path = NULL) {
  prompt <- is.null(path)
  if (prompt) {
    if (!interactive()) {
      rlang::abort(
        "Supply `path` to `set_output_dir()` in non-interactive sessions."
      )
    }
    path <- readline("Output directory for file-based runs (Enter to cancel): ")
    if (!nzchar(path)) {
      return(invisible(NULL))
    }
  }

  checkmate::assert_string(path, min.chars = 1L, pattern = "^[^\r\n]+$")
  path <- fs::path_abs(fs::path_expand(path))
  if (fs::file_exists(path) && !fs::dir_exists(path)) {
    rlang::abort("`path` must point to a directory.")
  }
  if (prompt) {
    answer <- readline(sprintf("Save '%s' for future runs? [y/N] ", path))
    if (!tolower(answer) %in% c("y", "yes")) {
      return(invisible(NULL))
    }
  }

  config <- output_dir_config_path()
  fs::dir_create(fs::path_dir(config))
  temporary <- tempfile("output-dir-", tmpdir = fs::path_dir(config))
  on.exit(unlink(temporary), add = TRUE)
  writeLines(enc2utf8(path), temporary, useBytes = TRUE)
  if (!file.rename(temporary, config)) {
    rlang::abort("Could not save the output directory preference.")
  }
  invisible(path)
}

output_dir_config_path <- function() {
  fs::path(tools::R_user_dir("immundapter", "config"), "output-dir")
}

read_saved_output_dir <- function() {
  config <- output_dir_config_path()
  if (fs::is_file(config)) {
    path <- readLines(config, n = 1L, encoding = "UTF-8", warn = FALSE)
    if (length(path) && fs::is_absolute_path(path)) {
      return(path)
    }
  }
  NULL
}
