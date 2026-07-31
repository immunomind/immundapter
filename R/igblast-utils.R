igblast_runtime_info <- function() {
  if (!suppressMessages(requireNamespace("igblastr", quietly = TRUE))) {
    return(list(
      package_installed = FALSE,
      igblast_installed = FALSE
    ))
  }

  installed <- isTRUE(getExportedValue("igblastr", "has_igblast")())
  runtime <- list(
    package_installed = TRUE,
    igblastr_version = as.character(utils::packageVersion("igblastr")),
    igblast_installed = installed
  )
  if (!installed) {
    return(runtime)
  }

  info <- getExportedValue("igblastr", "igblast_info")()
  runtime$igblast_version <- as.character(info$igblastn_version)
  runtime$igblast_root <- as.character(info$igblast_root)
  runtime$germline_db <- tryCatch(
    suppressMessages(getExportedValue("igblastr", "use_germline_db")()),
    error = function(error) NULL
  )
  runtime$c_region_db <- tryCatch(
    suppressMessages(getExportedValue("igblastr", "use_c_region_db")()),
    error = function(error) ""
  )
  runtime
}


igblast_prepare_input <- function(idata, sequence_col) {
  chain_col <- immundata::imd_schema("chain")
  required <- c(chain_col, sequence_col)
  missing <- setdiff(required, colnames(idata$annotations))
  if (length(missing) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST input is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  chains <- idata$annotations |>
    dplyr::select(dplyr::all_of(required)) |>
    dplyr::rename(sequence = dplyr::all_of(sequence_col)) |>
    dplyr::filter(!is.na(sequence)) |>
    dplyr::filter(sequence != "") |>
    dplyr::mutate(sequence = toupper(sequence)) |>
    dplyr::distinct()

  queries <- chains |>
    dplyr::select(sequence) |>
    dplyr::distinct() |>
    dplyr::arrange(sequence) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(queries) == 0L) {
    rlang::abort(sprintf(
      "No nucleotide sequences are available in `%s`.",
      sequence_col
    ))
  }
  if (
    !is.character(queries$sequence) ||
      any(!grepl("^[ACGTRYSWKMBDHVN]+$", queries$sequence))
  ) {
    rlang::abort(sprintf(
      "`%s` must contain ungapped IUPAC nucleotide sequences.",
      sequence_col
    ))
  }

  queries$igblast_query_id <- sprintf(
    "igblast_query_%09d",
    seq_len(nrow(queries))
  )
  queries <- queries[c("igblast_query_id", "sequence")]

  lookup <- duckplyr::as_duckdb_tibble(queries)
  mapping <- chains |>
    dplyr::inner_join(lookup, by = "sequence") |>
    dplyr::select(
      dplyr::all_of(c(chain_col, "igblast_query_id"))
    ) |>
    dplyr::distinct()

  list(queries = queries, mapping = mapping)
}


igblast_call <- function(arguments) {
  invisible(do.call(
    getExportedValue("igblastr", "igblastn"),
    arguments
  ))
}
