######################################
### Utility functions for TCRdist3 ###
######################################


tcrdist3_chain_columns <- function(chain) {
  suffix <- switch(chain,
    alpha = "a",
    beta = "b",
    gamma = "g",
    delta = "d"
  )
  c(
    cdr3 = paste0("cdr3_", suffix, "_aa"),
    v_gene = paste0("v_", suffix, "_gene"),
    j_gene = paste0("j_", suffix, "_gene")
  )
}

tcrdist3_chain_locus <- function(chain) {
  switch(chain,
    alpha = "TRA",
    beta = "TRB",
    gamma = "TRG",
    delta = "TRD"
  )
}

tcrdist3_rename_chain_columns <- function(table,
                                          chain,
                                          config,
                                          source_suffix = "") {
  target <- tcrdist3_chain_columns(chain)
  source <- paste0(c(
    cdr3 = config$cdr3_cols[[chain]],
    v_gene = config$v_gene_cols[[chain]],
    j_gene = config$j_gene_cols[[chain]]
  ), source_suffix)
  names(source) <- names(target)
  missing <- setdiff(unname(source), colnames(table))
  if (length(missing) > 0L) {
    rlang::abort(sprintf(
      "tcrdist3 input is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  for (column in unname(source)) {
    table <- table |>
      dplyr::filter(!is.na(!!rlang::sym(column))) |>
      dplyr::filter(!!rlang::sym(column) != "")
  }

  table |>
    dplyr::rename(
      !!target[["cdr3"]] := dplyr::all_of(source[["cdr3"]]),
      !!target[["v_gene"]] := dplyr::all_of(source[["v_gene"]]),
      !!target[["j_gene"]] := dplyr::all_of(source[["j_gene"]])
    )
}

tcrdist3_schema_chain_loci <- function(chain) {
  trimws(strsplit(chain, "|", fixed = TRUE)[[1L]])
}

tcrdist3_prepare_immundata <- function(input, config) {
  receptor_col <- immundata::imd_schema("receptor")
  locus_col <- immundata::imd_schema("locus")
  schema_chains <- immundata::imd_receptor_chains(input$schema_receptor)
  prepared <- input$receptors
  selected_columns <- receptor_col

  if (length(schema_chains) < 2L) {
    if (length(config$chains) != 1L) {
      rlang::abort(sprintf(
        paste(
          "The `ImmunData` receptor schema is not paired, but tcrdist3",
          "is configured for chains: %s."
        ),
        paste(config$chains, collapse = ", ")
      ))
    }

    chain <- config$chains[[1L]]
    if (
      length(schema_chains) == 1L &&
        !tcrdist3_chain_locus(chain) %in%
          tcrdist3_schema_chain_loci(schema_chains[[1L]])
    ) {
      rlang::abort(sprintf(
        paste(
          "tcrdist3 chain `%s` is incompatible with the",
          "`ImmunData` receptor schema chain `%s`."
        ),
        chain,
        schema_chains[[1L]]
      ))
    }

    prepared <- tcrdist3_rename_chain_columns(prepared, chain, config)
    selected_columns <- c(
      selected_columns,
      unname(tcrdist3_chain_columns(chain))
    )
  } else {
    for (chain in config$chains) {
      chain_locus <- tcrdist3_chain_locus(chain)
      schema_position <- which(vapply(
        schema_chains,
        function(schema_chain) {
          chain_locus %in% tcrdist3_schema_chain_loci(schema_chain)
        },
        logical(1L)
      ))
      if (length(schema_position) != 1L) {
        rlang::abort(sprintf(
          paste(
            "tcrdist3 chain `%s` is incompatible with the",
            "`ImmunData` paired receptor schema: %s."
          ),
          chain,
          paste(schema_chains, collapse = ", ")
        ))
      }

      side <- c("x", "y")[[schema_position]]
      source_suffix <- paste0(".", side)
      locus_source <- paste0(locus_col, source_suffix)
      if (!locus_source %in% colnames(prepared)) {
        rlang::abort(sprintf(
          "The `ImmunData` receptors are missing column `%s`.",
          locus_source
        ))
      }

      prepared <- prepared |>
        dplyr::filter(
          !!rlang::sym(locus_source) == !!chain_locus
        )
      prepared <- tcrdist3_rename_chain_columns(
        prepared,
        chain,
        config,
        source_suffix = source_suffix
      )
      selected_columns <- c(
        selected_columns,
        unname(tcrdist3_chain_columns(chain))
      )
    }
  }

  prepared |>
    dplyr::select(dplyr::all_of(selected_columns)) |>
    dplyr::rename(
      input_id = dplyr::all_of(receptor_col)
    )
}

tcrdist3_read_input_file <- function(path) {
  extension <- tolower(fs::path_ext(path))
  options <- if (extension %in% c("tsv", "txt")) {
    list(delim = "\t")
  } else {
    list()
  }

  duckplyr::read_csv_duckdb(path, options = options)
}

tcrdist3_prepare_file <- function(input, config) {
  table <- tcrdist3_read_input_file(input)
  source_columns <- unique(c(
    unname(config$cdr3_cols),
    unname(config$v_gene_cols),
    unname(config$j_gene_cols)
  ))
  missing <- setdiff(source_columns, colnames(table))
  if (length(missing) > 0L) {
    rlang::abort(sprintf(
      "tcrdist3 input is missing required column(s): %s.",
      paste(missing, collapse = ", ")
    ))
  }

  id_column <- if ("imd_receptor_id" %in% colnames(table)) {
    "imd_receptor_id"
  } else if ("cell_id" %in% colnames(table)) {
    "cell_id"
  } else {
    NULL
  }

  selected <- c(id_column, source_columns)
  table <- table |>
    dplyr::select(dplyr::all_of(selected))
  for (chain in config$chains) {
    table <- tcrdist3_rename_chain_columns(table, chain, config)
  }

  if (is.null(id_column)) {
    table <- table |>
      dplyr::mutate(input_id = dplyr::row_number())
  } else {
    table <- table |>
      dplyr::rename(input_id = dplyr::all_of(id_column))
  }

  table |>
    dplyr::distinct()
}

tcrdist3_prepare_input <- function(input, config) {
  if (inherits(input, "ImmunData")) {
    return(tcrdist3_prepare_immundata(input, config))
  }
  tcrdist3_prepare_file(input, config)
}
