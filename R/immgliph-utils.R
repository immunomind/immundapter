immgliph_runtime_info <- function() {
  installed <- requireNamespace("immGLIPH", quietly = TRUE) &&
    exists("runGLIPH", envir = asNamespace("immGLIPH"), inherits = FALSE)

  if (!installed) {
    return(list(installed = FALSE))
  }

  list(
    installed = TRUE,
    version = as.character(utils::packageVersion("immGLIPH")),
    path = system.file(package = "immGLIPH")
  )
}


immgliph_prepare_immundata <- function(input, config) {
  receptor_col <- immundata::imd_schema("receptor")
  locus_col <- immundata::imd_schema("locus")
  schema_chains <- immundata::imd_receptor_chains(input$schema_receptor)
  annotations <- input$annotations

  if (length(schema_chains) == 1L) {
    declared_loci <- trimws(strsplit(
      schema_chains[[1L]],
      "|",
      fixed = TRUE
    )[[1L]])
    if (!"TRB" %in% declared_loci) {
      rlang::abort(sprintf(
        "ImmGLIPH requires TRB data; the receptor schema contains `%s`.",
        schema_chains[[1L]]
      ))
    }
  } else {
    if (!locus_col %in% colnames(annotations)) {
      rlang::abort(sprintf(
        "ImmGLIPH input is missing the chain column `%s`.",
        locus_col
      ))
    }
    annotations <- annotations |>
      dplyr::filter(!!rlang::sym(locus_col) == "TRB")
  }

  input_columns <- c(
    CDR3b = config$cdr3_col,
    TRBV = config$v_gene_col,
    patient = config$patient_col,
    HLA = config$hla_col,
    counts = config$counts_col
  )
  selected_columns <- unique(c(receptor_col, unname(input_columns)))
  missing_columns <- setdiff(selected_columns, colnames(annotations))
  if (length(missing_columns) > 0L) {
    rlang::abort(sprintf(
      "ImmGLIPH input is missing required column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }

  mapping <- annotations |>
    dplyr::select(dplyr::all_of(selected_columns)) |>
    dplyr::filter(!is.na(!!rlang::sym(config$cdr3_col))) |>
    dplyr::filter(!!rlang::sym(config$cdr3_col) != "") |>
    dplyr::distinct()

  list(mapping = mapping, input_columns = input_columns)
}


immgliph_build_groups <- function(result, mapping, cdr3_col) {
  clusters <- result$cluster_list
  if (is.null(clusters)) {
    rlang::abort("ImmGLIPH output is missing `cluster_list`.")
  }

  if (length(clusters) == 0L) {
    membership <- data.frame(
      CDR3b = character(),
      immgliph_group_id = character()
    )
  } else {
    group_ids <- names(clusters)
    if (is.null(group_ids) || any(!nzchar(group_ids))) {
      group_ids <- sprintf("immgliph_group_%06d", seq_along(clusters))
    }

    membership <- do.call(rbind, Map(
      function(cluster, group_id) {
        if (is.null(cluster[["CDR3b"]])) {
          rlang::abort("An ImmGLIPH group is missing its `CDR3b` column.")
        }
        data.frame(
          CDR3b = as.character(cluster[["CDR3b"]]),
          immgliph_group_id = group_id,
          stringsAsFactors = FALSE
        )
      },
      clusters,
      group_ids
    ))
  }

  receptor_col <- immundata::imd_schema("receptor")
  membership <- duckplyr::as_duckdb_tibble(membership)

  mapping |>
    dplyr::select(dplyr::all_of(c(receptor_col, cdr3_col))) |>
    dplyr::rename(CDR3b = dplyr::all_of(cdr3_col)) |>
    dplyr::inner_join(membership, by = "CDR3b") |>
    dplyr::rename(dplyr::all_of(
      stats::setNames("CDR3b", cdr3_col)
    )) |>
    dplyr::select(dplyr::all_of(c(
      receptor_col,
      cdr3_col,
      "immgliph_group_id"
    ))) |>
    dplyr::distinct()
}


immgliph_call <- function(arguments) {
  do.call(getExportedValue("immGLIPH", "runGLIPH"), arguments)
}
