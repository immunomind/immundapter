raxml_ng_executable <- function() {
  Sys.which("raxml-ng")
}


raxml_ng_execute <- function(executable, args, stdout, stderr) {
  suppressWarnings(system2(
    executable,
    args,
    stdout = stdout,
    stderr = stderr
  ))
}


raxml_ng_command_args <- function(config,
                                  alignment_path,
                                  prefix,
                                  model_path = NULL) {
  checkmate::assert_r6(config, classes = "RaxmlNgToolConfig")
  checkmate::assert_string(alignment_path, min.chars = 1L)
  checkmate::assert_string(prefix, min.chars = 1L)
  checkmate::assert_string(model_path, min.chars = 1L, null.ok = TRUE)

  partitioned <- !is.null(config$partition_col)
  if (partitioned && is.null(model_path)) {
    rlang::abort("`model_path` is required for partitioned RAxML-NG input.")
  }
  if (!partitioned && !is.null(model_path)) {
    rlang::abort("`model_path` requires a configured `partition_col`.")
  }

  args <- c(
    paste0("--", config$analysis),
    "--msa", shQuote(alignment_path),
    "--model", shQuote(if (partitioned) model_path else config$model),
    "--prefix", shQuote(prefix)
  )
  if (!is.null(config$seed)) {
    args <- c(args, "--seed", as.character(config$seed))
  }
  if (!is.null(config$threads)) {
    args <- c(args, "--threads", as.character(config$threads))
  }
  if (partitioned && !is.null(config$brlen)) {
    args <- c(args, "--brlen", config$brlen)
  }
  args
}


raxml_ng_format_command <- function(executable, args) {
  paste(c(shQuote(executable), args), collapse = " ")
}


raxml_ng_parse_partitions <- function(specification, alignment_length) {
  checkmate::assert_string(specification, min.chars = 1L)
  checkmate::assert_count(alignment_length, positive = TRUE)

  if (grepl("(^|;)\\s*(;|$)", specification)) {
    rlang::abort("Partition specifications cannot contain empty entries.")
  }
  entries <- trimws(strsplit(specification, ";", fixed = TRUE)[[1L]])
  if (length(entries) == 0L || any(!nzchar(entries))) {
    rlang::abort("Partition specifications cannot contain empty entries.")
  }

  parsed <- lapply(entries, function(entry) {
    match <- regexec(
      "^([^=]+?)\\s*=\\s*([0-9]+)\\s*-\\s*([0-9]+)$",
      entry,
      perl = TRUE
    )
    fields <- regmatches(entry, match)[[1L]]
    if (length(fields) == 0L) {
      rlang::abort(sprintf(
        "Invalid partition entry `%s`; expected `name=start-end`.",
        entry
      ))
    }

    data.frame(
      name = trimws(fields[[2L]]),
      start = as.numeric(fields[[3L]]),
      end = as.numeric(fields[[4L]]),
      stringsAsFactors = FALSE
    )
  })
  partitions <- do.call(rbind, parsed)

  if (any(!nzchar(partitions$name)) ||
      any(grepl("[=;,\r\n]", partitions$name))) {
    rlang::abort("Partition names must be non-empty and cannot contain `=`, `;`, or `,`.")
  }
  if (anyDuplicated(partitions$name)) {
    duplicated_names <- unique(partitions$name[duplicated(partitions$name)])
    rlang::abort(sprintf(
      "Partition names must be unique; duplicated: %s.",
      paste(duplicated_names, collapse = ", ")
    ))
  }
  if (any(!is.finite(partitions$start)) ||
      any(!is.finite(partitions$end)) ||
      any(partitions$start < 1) ||
      any(partitions$end < 1)) {
    rlang::abort("Partition coordinates must be positive one-based integers.")
  }
  if (any(partitions$start > partitions$end)) {
    rlang::abort("Every partition must have `start <= end`.")
  }
  if (any(partitions$end > alignment_length)) {
    rlang::abort(sprintf(
      "Partition coordinates exceed the family alignment length (%d).",
      alignment_length
    ))
  }

  partitions <- partitions[order(partitions$start, partitions$end), , drop = FALSE]
  rownames(partitions) <- NULL
  if (nrow(partitions) > 1L &&
      any(partitions$start[-1L] <= partitions$end[-nrow(partitions)])) {
    rlang::abort("Partition coordinates cannot overlap.")
  }
  complete <- partitions$start[[1L]] == 1 &&
    partitions$end[[nrow(partitions)]] == alignment_length &&
    (nrow(partitions) == 1L || all(
      partitions$start[-1L] == partitions$end[-nrow(partitions)] + 1
    ))
  if (!complete) {
    rlang::abort(sprintf(
      "Partitions must cover every alignment position from 1 to %d.",
      alignment_length
    ))
  }

  partitions
}


raxml_ng_write_model_file <- function(partitions, model, path) {
  checkmate::assert_data_frame(
    partitions,
    min.rows = 1L,
    col.names = "named"
  )
  checkmate::assert_names(
    names(partitions),
    must.include = c("name", "start", "end")
  )
  checkmate::assert_string(model, min.chars = 1L)
  checkmate::assert_string(path, min.chars = 1L)

  lines <- sprintf(
    "%s, %s = %.0f-%.0f",
    model,
    partitions$name,
    partitions$start,
    partitions$end
  )
  writeLines(lines, path, useBytes = TRUE)
  invisible(path)
}


raxml_ng_prepare_immundata <- function(input, config) {
  locus_col <- immundata::imd_schema("locus")
  schema_chains <- immundata::imd_receptor_chains(input$schema_receptor)
  annotations <- input$annotations

  if (length(schema_chains) == 1L) {
    declared_loci <- trimws(strsplit(schema_chains[[1L]], "|", fixed = TRUE)[[1L]])

    if (!config$locus %in% declared_loci) {
      rlang::abort(sprintf(
        paste(
          "RAxML-NG locus `%s` is incompatible with the",
          "`ImmunData` receptor schema chain `%s`."
        ),
        config$locus,
        schema_chains[[1L]]
      ))
    }
  } else if (length(schema_chains) == 2L) {
    matching_chains <- which(vapply(
      schema_chains,
      function(chain) {
        config$locus %in% trimws(strsplit(chain, "|", fixed = TRUE)[[1L]])
      },
      logical(1L)
    ))
    if (length(matching_chains) != 1L) {
      rlang::abort(sprintf(
        paste(
          "RAxML-NG locus `%s` is incompatible with the",
          "`ImmunData` paired receptor schema: %s."
        ),
        config$locus,
        paste(schema_chains, collapse = ", ")
      ))
    }
    if (!locus_col %in% colnames(annotations)) {
      rlang::abort(sprintf(
        "The `ImmunData` annotations are missing column `%s`.",
        locus_col
      ))
    }
    annotations <- annotations |>
      dplyr::filter(!!rlang::sym(locus_col) == config$locus)
  } else if (locus_col %in% colnames(annotations)) {
    annotations <- annotations |>
      dplyr::filter(!!rlang::sym(locus_col) == config$locus)
  }

  source <- c(
    family_id = config$family_col,
    tip_id = config$tip_col,
    aligned_sequence = config$alignment_col
  )
  if (!is.null(config$partition_col)) {
    source <- c(source, partition_spec = config$partition_col)
  }
  missing_columns <- setdiff(unname(source), colnames(annotations))
  if (length(missing_columns) > 0L) {
    rlang::abort(sprintf(
      "`input$annotations` is missing RAxML-NG column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }

  annotations |>
    dplyr::select(dplyr::all_of(unname(source))) |>
    dplyr::rename(dplyr::all_of(source)) |>
    dplyr::distinct()
}


raxml_ng_prepare_additional_tips <- function(additional_tips,
                                             observed,
                                             config) {
  if (is.null(additional_tips)) {
    return(NULL)
  }
  if (!is.data.frame(additional_tips)) {
    rlang::abort("`additional_tips` must be a data frame or duckplyr table.")
  }

  required <- c("family_id", "tip_id", "aligned_sequence")
  missing_columns <- setdiff(required, colnames(additional_tips))
  if (length(missing_columns) > 0L) {
    rlang::abort(sprintf(
      "`additional_tips` is missing required column(s): %s.",
      paste(missing_columns, collapse = ", ")
    ))
  }

  columns <- c(
    required,
    if ("partition_spec" %in% colnames(additional_tips)) "partition_spec"
  )
  additional_tips <- additional_tips |>
    dplyr::select(dplyr::all_of(columns)) |>
    duckplyr::as_duckdb_tibble()

  observed_families <- observed |>
    dplyr::select(dplyr::all_of("family_id")) |>
    dplyr::distinct() |>
    dplyr::collect()
  observed_families <- as.character(observed_families$family_id)
  additional_families <- additional_tips |>
    dplyr::select(dplyr::all_of("family_id")) |>
    dplyr::distinct() |>
    dplyr::collect()
  additional_families <- as.character(additional_families$family_id)
  valid_family <- !is.na(additional_families) &
    nzchar(trimws(additional_families))
  unknown_families <- setdiff(
    additional_families[valid_family],
    observed_families
  )
  if (length(unknown_families) > 0L) {
    rlang::abort(sprintf(
      "`additional_tips` contains unknown family/families: %s.",
      paste(unknown_families, collapse = ", ")
    ))
  }

  if (!is.null(config$partition_col) &&
      !"partition_spec" %in% colnames(additional_tips)) {
    family_partitions <- observed |>
      dplyr::select(dplyr::all_of(c("family_id", "partition_spec"))) |>
      dplyr::distinct()
    additional_tips <- additional_tips |>
      dplyr::left_join(family_partitions, by = "family_id")
  }
  if (is.null(config$partition_col) &&
      "partition_spec" %in% colnames(additional_tips)) {
    additional_tips <- additional_tips |>
      dplyr::select(-dplyr::all_of("partition_spec"))
  }

  additional_tips
}


raxml_ng_collect_input <- function(input, config, additional_tips = NULL) {
  observed <- raxml_ng_prepare_immundata(input, config)
  additional_tips <- raxml_ng_prepare_additional_tips(
    additional_tips,
    observed,
    config
  )

  if (!is.null(additional_tips)) {
    observed <- observed |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    additional_tips <- additional_tips |>
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    observed <- dplyr::bind_rows(observed, additional_tips)
  }

  prepared <- observed |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(prepared) == 0L) {
    rlang::abort("No receptors are available for the configured RAxML-NG locus.")
  }

  for (column in c(
    "family_id",
    "tip_id",
    "aligned_sequence",
    if (!is.null(config$partition_col)) "partition_spec"
  )) {
    values <- as.character(prepared[[column]])
    invalid <- is.na(prepared[[column]]) | !nzchar(trimws(values))
    if (any(invalid)) {
      rlang::abort(sprintf(
        "RAxML-NG input contains %d missing or empty `%s` value(s).",
        sum(invalid),
        column
      ))
    }
    if (any(grepl("[\t\r\n]", values))) {
      rlang::abort(sprintf(
        "RAxML-NG input `%s` values cannot contain tabs or line breaks.",
        column
      ))
    }
    prepared[[column]] <- values
  }
  if (any(grepl("[[:space:]]", prepared$aligned_sequence))) {
    rlang::abort(
      "RAxML-NG aligned sequences cannot contain whitespace characters."
    )
  }
  prepared$aligned_sequence <- toupper(prepared$aligned_sequence)
  invalid_sequence <- grepl(
    "[^ACGTURYSWKMBDHVNXO?.-]",
    prepared$aligned_sequence
  )
  if (any(invalid_sequence)) {
    invalid_chars <- unique(strsplit(
      gsub(
        "[ACGTURYSWKMBDHVNXO?.-]",
        "",
        paste(prepared$aligned_sequence[invalid_sequence], collapse = "")
      ),
      ""
    )[[1L]])
    rlang::abort(sprintf(
      "RAxML-NG aligned sequences contain unsupported character(s): %s.",
      paste(invalid_chars, collapse = ", ")
    ))
  }
  prepared <- unique(prepared)

  if (!is.null(config$partition_col)) {
    specifications <- split(prepared$partition_spec, prepared$family_id)
    conflicting <- vapply(
      specifications,
      function(values) length(unique(values)) != 1L,
      logical(1L)
    )
    if (any(conflicting)) {
      rlang::abort(sprintf(
        "Conflicting partition specifications in family/families: %s.",
        paste(names(conflicting)[conflicting], collapse = ", ")
      ))
    }
  }

  duplicate_tip <- duplicated(prepared[c("family_id", "tip_id")])
  if (any(duplicate_tip)) {
    example <- prepared[which(duplicate_tip)[[1L]], , drop = FALSE]
    rlang::abort(sprintf(
      paste(
        "RAxML-NG tip `%s` has conflicting aligned sequences",
        "in family `%s`."
      ),
      example$tip_id,
      example$family_id
    ))
  }

  lengths_per_family <- split(
    nchar(prepared$aligned_sequence),
    prepared$family_id
  )
  aligned <- vapply(lengths_per_family, function(lengths) {
    length(unique(lengths)) == 1L
  }, logical(1L))
  if (any(!aligned)) {
    rlang::abort(sprintf(
      "RAxML-NG sequences are not aligned to one length in family/families: %s.",
      paste(names(aligned)[!aligned], collapse = ", ")
    ))
  }

  if (!is.null(config$partition_col)) {
    for (family_id in names(lengths_per_family)) {
      raxml_ng_parse_partitions(
        specifications[[family_id]][[1L]],
        unique(lengths_per_family[[family_id]])[[1L]]
      )
    }
  }

  prepared[order(prepared$family_id, prepared$tip_id), , drop = FALSE]
}
