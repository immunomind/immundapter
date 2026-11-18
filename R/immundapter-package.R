#' immundapter: Run external immune-repertoire tools from ImmunData
#'
#' immundapter connects immune-repertoire data in R to external analysis tools
#' and reference databases. Each tool prepares the required input, runs the
#' external software, records the generated files, and provides results as
#' tables, file paths, or updated [immundata::ImmunData] datasets.
#'
#' @section Analysis-tool workflow:
#' Analysis tools use the same sequence:
#' ```
#' config <- conf_tcrdist3()
#' tool <- tool_tcrdist3(config)
#' run <- tool$run(input, output_dir = "analysis/tcrdist3")
#' results <- run$results()
#' ```
#'
#' The accepted `input` is different for each tool. See [tool_tcrdist3()],
#' [tool_igblast()], [immgliph], or [raxml-ng] for exact input columns and
#' result fields. The `$run()` records where the files were saved. Use
#' `$results()` to read the biological output.
#'
#' Use [immundapter::set_output_dir()] to save a default for file-based
#' runs. An explicit `output_dir` overrides the saved preference for that run.
#'
#' Large result tables remain on disk until needed. Use `dplyr::collect()` to
#' read a table into R.
#'
#' VDJdb is a reference database, so it uses a shorter workflow.
#' See [data_vdjdb()] for downloading, loading, and matching VDJdb records.
#'
#' @keywords internal
#' @importFrom R6 R6Class
#' @importFrom checkmate assert_choice assert_class assert_data_frame assert_file_exists assert_int assert_list assert_r6 assert_string
#' @importFrom cli cli_inform
#' @importFrom dplyr all_of collect distinct filter inner_join rename select
#' @importFrom duckplyr as_duckdb_tibble compute_parquet read_parquet_duckdb
#' @importFrom fs dir_create path path_abs
#' @importFrom immundata annotate_receptors imd_schema
#' @importFrom jsonlite write_json
#' @importFrom rlang abort sym
#' @importFrom withr with_seed
"_PACKAGE"

utils::globalVariables(c(
  "dd", "v_call", "germline_alignment", ".position", ".germline_base",
  "sequence_alignment", "sequence_id", ".query_base",
  ".trimmed_sequence", "v_germline_start", "pattern",
  ".normalized_sequence", "replacement", "imgt_v_end",
  "v_germline_end", "j_germline_start", "j_germline_end"
))
