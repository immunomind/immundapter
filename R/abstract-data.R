#' Base class for immune-repertoire data sources
#'
#' `AbstractData` contains the code shared by reference-data sources.
#'
#' Data are stored in `~/immundapter_data/data` by default. Set
#' `IMMUNDAPTER_DATA` to change the root directory. A concrete data-source
#' constructor may accept another output directory.
#'
#' @keywords internal
immundapter_data_root <- function() {
  root <- Sys.getenv("IMMUNDAPTER_DATA", unset = "")
  if (!nzchar(root)) {
    root <- fs::path_home("immundapter_data")
  }

  fs::path_abs(root)
}

immundapter_data_path <- function(...) {
  fs::path(immundapter_data_root(), ...)
}

AbstractData <- R6::R6Class(
  classname = "AbstractData",
  cloneable = FALSE,
  public = list(
    initialize = function(
      name,
      output_dir = immundapter_data_path("data")
    ) {
      checkmate::assert_string(name)
      checkmate::assert_string(output_dir)

      private$.name <- name
      private$.output_dir <- fs::path_abs(output_dir)
    },

    download = function(execute = FALSE, force = FALSE) {
      checkmate::assert_flag(execute)
      checkmate::assert_flag(force)
      private$download_impl(
        output_dir = private$.output_dir,
        execute = execute,
        force = force
      )
    },

    check = function() {
      private$check_impl(private$.output_dir)
    },

    load = function(output = c("duckplyr_df", "df")) {
      output <- match.arg(output)
      result <- private$load_impl(private$.output_dir)

      if (identical(output, "df")) {
        return(as.data.frame(dplyr::collect(result)))
      }
      result
    },

    annotate = function(idata, ...) {
      checkmate::assert_r6(idata, classes = "ImmunData")
      private$annotate_impl(idata, private$.output_dir, ...)
    },

    citation = function() {
      private$citation_impl()
    }
  ),
  active = list(
    name = function(value) {
      if (!missing(value)) {
        rlang::abort("`name` is read-only.")
      }
      private$.name
    },
    output_dir = function(value) {
      if (!missing(value)) {
        rlang::abort("`output_dir` is read-only.")
      }
      private$.output_dir
    }
  ),
  private = list(
    .name = NULL,
    .output_dir = NULL,

    download_impl = function(output_dir, execute, force) {
      rlang::abort("The data source must implement `download()`.")
    },
    check_impl = function(output_dir) {
      rlang::abort("The data source must implement `check()`.")
    },
    load_impl = function(output_dir) {
      rlang::abort("The data source must implement `load()`.")
    },
    annotate_impl = function(idata, output_dir, ...) {
      rlang::abort("The data source must implement `annotate()`.")
    },
    citation_impl = function() {
      rlang::abort("The data source must implement `citation()`.")
    }
  )
)
