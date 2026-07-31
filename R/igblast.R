#' Assign immune-receptor germline genes with IgBLAST
#'
#' `tool_igblast()` runs `igblastr::igblastn()` directly on chain nucleotide
#' sequences from an [immundata::ImmunData] object. The germline database
#' selected with `igblastr::use_germline_db()` determines both the organism and
#' whether Ig or TCR sequences are analyzed.
#'
#' Output is always AIRR, domains are always IMGT, and igblastr automatically
#' resolves germline, optional C-region, internal, and auxiliary data. Exact
#' duplicate sequences are submitted once and related back to `imd_chain_id`
#' through a persisted query map.
#'
#' @section Result status:
#' Results expose the complete AIRR output and query mapping as lazy duckplyr
#' tables. `$annotate()` adds selected AIRR columns to chain annotations,
#' defaulting to `sequence_alignment` and `germline_alignment`.
#'
#' @param sequence_col Column in `idata$annotations` containing ungapped
#'   nucleotide query sequences.
#' @param num_threads Optional number of IgBLAST threads. `NULL` uses IgBLAST's
#'   default.
#' @param run_args Additional named arguments passed through `...` in
#'   `igblastr::igblastn()`. Adapter-owned arguments cannot be overridden.
#' @param config An [IgblastToolConfig] created by [conf_igblast()].
#'
#' @return [conf_igblast()] returns a read-only [IgblastToolConfig].
#'   [tool_igblast()] returns a reusable IgBLAST adapter.
#'
#' @examples
#' config <- conf_igblast(sequence_col = "sequence")
#' tool <- tool_igblast(config)
#' tool$help()
#'
#' \dontrun{
#' pak::pak("bioc::igblastr")
#' igblastr::install_igblast()
#' igblastr::install_IMGT_germline_db(
#'   "<release>",
#'   organism = "Homo sapiens",
#'   tcr.db = TRUE
#' )
#' igblastr::use_germline_db("<db-name>")
#'
#' run <- tool$run(idata, output_dir = "analysis/igblast")
#' results <- run$results()
#' dplyr::collect(results$airr)
#' dplyr::collect(results$mapping)
#' aligned_idata <- results$annotate(idata)
#' }
#'
#' @seealso
#'   \url{https://bioconductor.org/packages/release/bioc/vignettes/igblastr/inst/doc/igblastr_overview.html}
#' @export
conf_igblast <- function(sequence_col = "sequence",
                         num_threads = NULL,
                         run_args = list()) {
  IgblastToolConfig$new(list(
    sequence_col = sequence_col,
    num_threads = num_threads,
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
#' @field sequence_col Column containing ungapped nucleotide query sequences.
#' @field num_threads Optional number of IgBLAST threads.
#' @field run_args Additional named IgBLAST command-line arguments.
#'
#' @keywords internal
IgblastToolConfig <- R6::R6Class(
  classname = "IgblastToolConfig",
  inherit = AbstractConfig,
  cloneable = FALSE,
  public = list(
    sequence_col = NULL,
    num_threads = NULL,
    run_args = NULL,
    #' @description
    #' Create a validated IgBLAST configuration. Users should call
    #' [conf_igblast()] instead.
    #' @param values A named list of IgBLAST configuration values.
    initialize = function(values) {
      checkmate::assert_list(values, names = "unique")
      checkmate::assert_string(values$sequence_col)
      checkmate::assert_int(values$num_threads, lower = 1L, null.ok = TRUE)
      checkmate::assert_list(values$run_args)

      argument_names <- names(values$run_args)
      if (
        length(values$run_args) > 0L &&
          (
            is.null(argument_names) ||
              any(!nzchar(argument_names)) ||
              anyDuplicated(argument_names)
          )
      ) {
        rlang::abort("IgBLAST `run_args` must have unique, non-empty names.")
      }

      fixed <- c(
        "query", "outfmt", "germline_db_V", "germline_db_D",
        "germline_db_J", "organism", "c_region_db",
        "custom_internal_data", "auxiliary_data", "domain_system",
        "ig_seqtype", "out", "parse.out", "show.in.browser",
        "show.command.only", "num_threads"
      )
      overridden <- intersect(argument_names, fixed)
      if (length(overridden) > 0L) {
        rlang::abort(sprintf(
          "IgBLAST `run_args` cannot override: %s.",
          paste(overridden, collapse = ", ")
        ))
      }

      super$initialize(values)
    }
  )
)


#' IgBLAST adapter
#'
#' @keywords internal
IgblastTool <- R6::R6Class(
  classname = "IgblastTool",
  inherit = AbstractTool,
  cloneable = FALSE,
  public = list(
    #' @description
    #' Create an IgBLAST adapter. Users should call [tool_igblast()] instead.
    #' @param config An [IgblastToolConfig] returned by [conf_igblast()].
    initialize = function(config = conf_igblast()) {
      checkmate::assert_r6(config, classes = "IgblastToolConfig")
      super$initialize("igblast", config, IgblastRun)
    }
  ),
  private = list(
    help_impl = function() {
      cli::cli_inform(c(
        "i" = "Set up IgBLAST through {.pkg igblastr}:",
        ">" = '{.code pak::pak("bioc::igblastr")}',
        ">" = "{.code igblastr::install_igblast()}",
        ">" = "{.code igblastr::list_germline_dbs()}",
        ">" = paste(
          "{.code igblastr::install_IMGT_germline_db(",
          '"<release>", organism = "Homo sapiens", tcr.db = TRUE)}'
        ),
        ">" = '{.code igblastr::use_germline_db("<db-name>")}',
        "i" = paste(
          "The selected database determines organism and Ig versus TCR.",
          "Selection must happen in the R session that calls `$run()`."
        ),
        "i" = paste(
          "A C-region database is optional. Install and select one with",
          "{.code install_IMGT_c_region_db()} and",
          "{.code use_c_region_db()} when constant-region calls are needed."
        )
      ))
    },
    check_impl = function() {
      runtime <- igblast_runtime_info()
      if (!runtime$package_installed) {
        cli::cli_inform(c(
          "x" = "{.pkg igblastr} is not installed.",
          "i" = "Run {.code tool$help()} for setup instructions."
        ))
        return(FALSE)
      }
      if (!runtime$igblast_installed) {
        cli::cli_inform(c(
          "v" = "{.pkg igblastr} {runtime$igblastr_version} is installed.",
          "x" = "No working IgBLAST installation was found."
        ))
        return(FALSE)
      }
      if (is.null(runtime$germline_db)) {
        cli::cli_inform(c(
          "v" = "IgBLAST {runtime$igblast_version} is available.",
          "x" = "No germline V/D/J database is selected in this R session."
        ))
        return(FALSE)
      }

      cli::cli_inform(c(
        "v" = "{.pkg igblastr} {runtime$igblastr_version} is installed.",
        "v" = "IgBLAST {runtime$igblast_version} is available.",
        "v" = "Selected germline database: {runtime$germline_db}",
        "i" = if (nzchar(runtime$c_region_db)) {
          "Selected C-region database: {runtime$c_region_db}"
        } else {
          "No C-region database is selected."
        }
      ))
      TRUE
    },
    citation_impl = function() {
      rlang::abort("IgBLAST `citation()` is not implemented yet.")
    },
    validate_input_impl = function(input, config) {
      checkmate::assert_r6(input, classes = "ImmunData")
      input
    },
    run_impl = function(input, config, run_dir) {
      runtime <- igblast_runtime_info()
      if (!runtime$package_installed || !runtime$igblast_installed) {
        rlang::abort(
          "igblastr and a working IgBLAST installation are required."
        )
      }
      if (is.null(runtime$germline_db)) {
        rlang::abort(paste(
          "No germline V/D/J database is selected.",
          "Call `igblastr::use_germline_db()` in this R session."
        ))
      }

      prepared <- igblast_prepare_input(input, config$sequence_col)
      input_path <- fs::path(run_dir, "igblast-input.fasta")
      mapping_path <- fs::path(run_dir, "igblast-query-map.parquet")
      output_path <- fs::path(run_dir, "igblast-airr.tsv")
      config_path <- fs::path(run_dir, "igblast-config.json")

      writeLines(
        as.vector(rbind(
          paste0(">", prepared$queries$igblast_query_id),
          prepared$queries$sequence
        )),
        input_path,
        useBytes = TRUE
      )
      duckplyr::compute_parquet(prepared$mapping, mapping_path)

      fixed <- list(
        query = as.character(fs::path_abs(input_path)),
        outfmt = "AIRR",
        germline_db_V = "auto",
        germline_db_D = "auto",
        germline_db_J = "auto",
        organism = "auto",
        c_region_db = "auto",
        custom_internal_data = "auto",
        auxiliary_data = "auto",
        domain_system = "imgt",
        ig_seqtype = "auto",
        out = as.character(fs::path_abs(output_path)),
        parse.out = FALSE,
        show.in.browser = FALSE,
        show.command.only = FALSE
      )
      arguments <- c(
        fixed,
        if (is.null(config$num_threads)) {
          list()
        } else {
          list(num_threads = config$num_threads)
        },
        config$run_args
      )

      jsonlite::write_json(
        list(
          sequence_col = config$sequence_col,
          num_threads = config$num_threads,
          run_args = config$run_args,
          fixed = fixed[setdiff(names(fixed), c("query", "out"))],
          runtime = runtime
        ),
        path = config_path,
        auto_unbox = TRUE,
        null = "null",
        pretty = TRUE
      )

      igblast_call(arguments)
      if (!fs::file_exists(output_path)) {
        rlang::abort("IgBLAST did not produce AIRR output.")
      }

      list(
        input = input_path,
        config = config_path,
        mapping = mapping_path,
        airr = output_path
      )
    }
  )
)
