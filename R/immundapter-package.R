#' immundapter: External tool adapters for ImmunData
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
