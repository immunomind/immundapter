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
  )
)
