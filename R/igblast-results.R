#' @noRd
IgblastResults <- R6::R6Class(
  classname = "IgblastResults",
  inherit = ToolResults,
  cloneable = FALSE,
  active = list(
    airr = function(value) {
      private$artifact(
        "airr",
        value,
        type = "csv",
        options = list(
          delim = "\t",
          header = TRUE,
          types = list(c(sequence_id = "VARCHAR"))
        )
      )
    },
    mapping = function(value) {
      private$artifact("mapping", value, type = "parquet")
    },
    raw_output = function(value) {
      private$artifact(
        "airr",
        value,
        type = "path",
        field = "raw_output"
      )
    },
    final = function(value) {
      private$artifact(
        "final",
        value,
        type = "csv",
        options = list(
          delim = "\t",
          header = TRUE,
          types = list(c(sequence_id = "VARCHAR"))
        ),
        missing_message = paste(
          "No final merged table is available.",
          "Run IgBLAST with named `fasta` and `annotations` input files."
        )
      )
    },
    final_output = function(value) {
      private$artifact(
        "final",
        value,
        type = "path",
        missing_message = paste(
          "No final merged table is available.",
          "Run IgBLAST with named `fasta` and `annotations` input files."
        ),
        field = "final_output"
      )
    }
  ),
  private = list(
    annotate_impl = function(
      idata,
      columns = c("sequence_alignment", "germline_alignment", "c_call"),
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
      if ("c_call" %in% columns && !"c_call" %in% colnames(airr)) {
        airr <- airr |>
          dplyr::mutate(c_call = NA_character_)
      }
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
