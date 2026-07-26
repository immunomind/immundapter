#' Configure or create an IgBLAST adapter
#'
#' `conf_igblast()` validates IgBLAST parameters. `tool_igblast()` creates the
#' adapter from that configuration.
#'
#' @param sequence_col Column containing nucleotide query sequences.
#' @param organism Query organism or `"auto"`.
#' @param ig_seqtype Receptor sequence type: `"Ig"`, `"TCR"`, or `"auto"`.
#' @param germline_db_V V-region database path or `"auto"`.
#' @param germline_db_D D-region database path or `"auto"`.
#' @param germline_db_J J-region database path or `"auto"`.
#' @param c_region_db Constant-region database path, `"auto"`, or `NULL`.
#' @param auxiliary_data Auxiliary-data file path, `"auto"`, or `NULL`.
#' @param num_threads Number of IgBLAST threads.
#' @param igblast_root Optional path to an existing IgBLAST installation.
#' @param run_args Additional named arguments passed to `igblastr::igblastn()`.
#' @param config An [IgblastToolConfig] created by `conf_igblast()`.
#'
#' @return `conf_igblast()` returns an IgBLAST configuration.
#'   `tool_igblast()` returns an IgBLAST adapter skeleton.
#'
#' @seealso \url{https://bioconductor.org/packages/release/bioc/manuals/igblastr/man/igblastr.pdf}
#' @export
conf_igblast <- function(sequence_col = "sequence",
                         organism = "auto",
                         ig_seqtype = "auto",
                         germline_db_V = "auto",
                         germline_db_D = "auto",
                         germline_db_J = "auto",
                         c_region_db = "auto",
                         auxiliary_data = "auto",
                         num_threads = 1L,
                         igblast_root = NULL,
                         run_args = list()) {
  IgblastToolConfig$new(list(
    sequence_col = sequence_col,
    organism = organism,
    ig_seqtype = ig_seqtype,
    germline_db_V = germline_db_V,
    germline_db_D = germline_db_D,
    germline_db_J = germline_db_J,
    c_region_db = c_region_db,
    auxiliary_data = auxiliary_data,
    num_threads = num_threads,
    igblast_root = igblast_root,
    run_args = run_args
  ))
}


#' @rdname conf_igblast
#' @export
tool_igblast <- function(config = conf_igblast()) {
  IgblastTool$new(config)
}


#' IgBLAST tool configuration
#'
#' @keywords internal
IgblastToolConfig <- R6::R6Class(
  classname = "IgblastToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    sequence_col = NULL,
    organism = NULL,
    ig_seqtype = NULL,
    germline_db_V = NULL,
    germline_db_D = NULL,
    germline_db_J = NULL,
    c_region_db = NULL,
    auxiliary_data = NULL,
    num_threads = NULL,
    igblast_root = NULL,
    run_args = NULL,
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_string(values$sequence_col)
      checkmate::assert_choice(
        values$organism,
        c("auto", "human", "mouse", "rat", "rabbit", "rhesus_monkey")
      )
      checkmate::assert_choice(values$ig_seqtype, c("auto", "Ig", "TCR"))
      checkmate::assert_string(values$germline_db_V)
      checkmate::assert_string(values$germline_db_D)
      checkmate::assert_string(values$germline_db_J)
      checkmate::assert_string(values$c_region_db, null.ok = TRUE)
      checkmate::assert_string(values$auxiliary_data, null.ok = TRUE)
      checkmate::assert_int(values$num_threads, lower = 1L)
      checkmate::assert_string(values$igblast_root, null.ok = TRUE)
      checkmate::assert_list(values$run_args, names = "unique")

      super$initialize(values)
    }
  )
)


#' IgBLAST adapter skeleton
#'
#' @keywords internal
IgblastTool <- R6::R6Class(
  classname = "IgblastTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    initialize = function(config = conf_igblast()) {
      checkmate::assert_r6(config, classes = "IgblastToolConfig")
      super$initialize("igblast", config)
    }
  ),
  private = list(
    help_impl = function() {
      install_package <- 'pak::pak("bioc::igblastr")'
      install_igblast <- "igblastr::install_igblast()"
      list_databases <- "igblastr::list_germline_dbs()"
      use_database <- 'igblastr::use_germline_db("<db-name>")'
      documentation <- paste0(
        "https://bioconductor.org/packages/release/bioc/",
        "vignettes/igblastr/inst/doc/igblastr_overview.html"
      )

      cli::cli_inform(c(
        "i" = "Install IgBLAST through {.pkg igblastr}:"
      ))
      div <- cli::cli_div(theme = list(
        ol = list("margin-left" = 2),
        ul = list("margin-left" = 2)
      ))
      steps <- cli::cli_ol(.close = FALSE)

      cli::cli_li("Install the R package on macOS, Linux, or Windows.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {install_package}}")
      cli::cli_end(commands)

      cli::cli_li("Install the platform-specific IgBLAST binary.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {install_igblast}}")
      cli::cli_end(commands)

      cli::cli_li("Choose a germline database.")
      commands <- cli::cli_ul(.close = FALSE)
      cli::cli_li("{.code {list_databases}}")
      cli::cli_li("{.code {use_database}}")
      cli::cli_end(commands)

      cli::cli_end(steps)
      cli::cli_end(div)

      cli::cli_inform(c("i" = "Documentation:"))
      div <- cli::cli_div(theme = list(
        ul = list("margin-left" = 2)
      ))
      links <- cli::cli_ul(.close = FALSE)
      cli::cli_li("igblastr overview: {.url {documentation}}")
      cli::cli_end(links)
      cli::cli_end(div)
    },
    check_impl = function() {
      package_installed <- requireNamespace("igblastr", quietly = TRUE)
      igblast_installed <- FALSE

      if (package_installed) {
        has_igblast <- getExportedValue("igblastr", "has_igblast")
        igblast_installed <- if (is.null(private$.config$igblast_root)) {
          isTRUE(has_igblast())
        } else {
          withr::with_envvar(
            c(IGBLAST_ROOT = private$.config$igblast_root),
            isTRUE(has_igblast())
          )
        }
      }

      if (igblast_installed) {
        version <- as.character(utils::packageVersion("igblastr"))
        path <- system.file(package = "igblastr")
        cli::cli_inform(c(
          "v" = "{.pkg igblastr} {version} is installed.",
          ">" = "R package: {.path {path}}",
          "v" = "A working IgBLAST installation is available."
        ))
      } else if (package_installed) {
        version <- as.character(utils::packageVersion("igblastr"))
        path <- system.file(package = "igblastr")
        cli::cli_inform(c(
          "v" = "{.pkg igblastr} {version} is installed.",
          ">" = "R package: {.path {path}}",
          "x" = "No working IgBLAST installation was found.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      } else {
        cli::cli_inform(c(
          "x" = "{.pkg igblastr} is not installed.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
      }

      package_installed && igblast_installed
    },
    citation_impl = function() {
      rlang::abort("IgBLAST `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      checkmate::assert_string(input, min.chars = 1L)
      checkmate::assert_file_exists(input)
      fs::path_abs(input)
    },
    run_impl = function(input, config, run_dir) {
      rlang::abort("IgBLAST `run()` is not implemented yet.")
    },
    results_impl = function(run) {
      rlang::abort("IgBLAST `results()` is not implemented yet.")
    },
    annotate_impl = function(idata, results, ...) {
      rlang::abort("IgBLAST `annotate()` is not implemented yet.")
    }
  )
)
