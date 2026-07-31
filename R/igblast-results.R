#' Parsed IgBLAST results
#'
#' `IgblastResults` preserves the complete AIRR output produced by IgBLAST.
#' `$annotate()` adds selected AIRR columns to chain-level `ImmunData`
#' annotations. By default it adds `sequence_alignment` and
#' `germline_alignment`; pass `columns` to select other AIRR fields.
#'
#' @section Annotation:
#' Use `results$annotate(idata, columns = ...)` to copy AIRR fields onto the
#' chains represented by the saved query mapping. Column names are preserved.
#' Additional arguments, including `conflicts`, are passed to
#' [immundata::annotate_chains()].
#'
#' @field airr Read-only lazy duckplyr table containing complete AIRR output.
#' @field mapping Read-only lazy duckplyr table mapping input `imd_chain_id`
#'   values to IgBLAST query identifiers.
#' @field raw_output Read-only path to the AIRR TSV file.
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

      airr <- self$airr
      missing <- setdiff(columns, colnames(airr))
      if (length(missing) > 0L) {
        rlang::abort(sprintf(
          "IgBLAST result column(s) not found: %s.",
          paste(missing, collapse = ", ")
        ))
      }

      chain_col <- immundata::imd_schema("chain")
      mapping <- self$mapping |>
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
