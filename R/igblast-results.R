#' Parsed IgBLAST results
#'
#' `IgblastResults` provides the complete AIRR output and the files created by
#' IgBLAST. `$airr`, `$mapping`, and `$final` stay on disk until you use
#' `dplyr::collect()`.
#'
#' @section Output depends on the input:
#' All runs provide `$airr`, `$mapping`, and `$raw_output`. Runs made from
#' `ImmunData` can use `results$annotate(idata, columns = ...)` to add AIRR
#' fields to the original chains. By default, it adds `sequence_alignment` and
#' `germline_alignment`.
#'
#' Runs made from a FASTA file cannot annotate `ImmunData` directly. Runs made
#' from a named FASTA-and-annotations pair also provide `$final` and
#' `$final_output`, containing the source table joined with AIRR results.
#'
#' Additional `$annotate()` arguments, including `conflicts`, are passed to
#' [immundata::annotate_chains()].
#'
#' @field airr Read-only AIRR result table.
#' @field mapping Read-only table showing which input chain or FASTA sequence
#'   produced each IgBLAST `sequence_id`.
#' @field raw_output Read-only file path to the AIRR TSV output.
#' @field final Read-only table containing the original annotation columns and
#'   AIRR results for paired FASTA and annotation input.
#' @field final_output Read-only file path to the final TSV table for paired
#'   input.
#'
#' @keywords internal
IgblastResults <- R6::R6Class(
  classname = "IgblastResults",
  inherit = ToolResults,
  cloneable = FALSE,
  active = list(
    airr = function(value) {
      if (!missing(value)) {
        rlang::abort("`airr` is read-only.")
      }
      path <- private$.run$artifacts$airr
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(
        path,
        options = list(delim = "\t", header = TRUE)
      )
    },
    mapping = function(value) {
      if (!missing(value)) {
        rlang::abort("`mapping` is read-only.")
      }
      path <- private$.run$artifacts$mapping
      checkmate::assert_file_exists(path)
      duckplyr::read_parquet_duckdb(path)
    },
    raw_output = function(value) {
      if (!missing(value)) {
        rlang::abort("`raw_output` is read-only.")
      }
      path <- private$.run$artifacts$airr
      checkmate::assert_file_exists(path)
      path
    },
    final = function(value) {
      if (!missing(value)) {
        rlang::abort("`final` is read-only.")
      }
      path <- private$.run$artifacts$final
      if (is.null(path)) {
        rlang::abort(paste(
          "No final merged table is available.",
          "Run IgBLAST with named `fasta` and `annotations` input files."
        ))
      }
      checkmate::assert_file_exists(path)
      duckplyr::read_csv_duckdb(
        path,
        options = list(delim = "\t", header = TRUE)
      )
    },
    final_output = function(value) {
      if (!missing(value)) {
        rlang::abort("`final_output` is read-only.")
      }
      path <- private$.run$artifacts$final
      if (is.null(path)) {
        rlang::abort(paste(
          "No final merged table is available.",
          "Run IgBLAST with named `fasta` and `annotations` input files."
        ))
      }
      checkmate::assert_file_exists(path)
      path
    }
  ),
  private = list(
    annotate_impl = function(
      idata,
      columns = c("sequence_alignment", "germline_alignment"),
      ...
    ) {
      checkmate::assert_character(
        columns,
        min.len = 1L,
        any.missing = FALSE,
        unique = TRUE
      )
      if ("sequence_id" %in% columns) {
        rlang::abort(
          "`sequence_id` identifies IgBLAST queries and cannot be annotated."
        )
      }

      chain_col <- immundata::imd_schema("chain")
      mapping <- self$mapping
      if (!chain_col %in% colnames(mapping)) {
        rlang::abort(paste(
          "IgBLAST results produced from FASTA input cannot annotate",
          "`ImmunData` because the input came from a FASTA file.",
          "Use `$airr` and `$mapping` to inspect or join the results."
        ))
      }

      airr <- self$airr
      missing <- setdiff(columns, colnames(airr))
      if (length(missing) > 0L) {
        rlang::abort(sprintf(
          "IgBLAST result column(s) not found: %s.",
          paste(missing, collapse = ", ")
        ))
      }

      mapping <- mapping |>
        dplyr::rename(
          sequence_id = igblast_query_id
        )
      annotations <- airr |>
        dplyr::select(
          dplyr::all_of(c("sequence_id", columns))
        ) |>
        dplyr::inner_join(mapping, by = "sequence_id") |>
        dplyr::select(
          dplyr::all_of(c(chain_col, columns))
        )

      immundata::annotate_chains(
        idata = idata,
        annotations = annotations,
        annot_col = chain_col,
        ...
      )
    }
  )
)
