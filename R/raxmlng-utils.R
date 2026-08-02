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


raxml_ng_collect_input <- function(input, config) {
  prepared <- raxml_ng_prepare_immundata(input, config) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(prepared) == 0L) {
    rlang::abort("No receptors are available for the configured RAxML-NG locus.")
  }

  for (column in c("family_id", "tip_id", "aligned_sequence")) {
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

  prepared[order(prepared$family_id, prepared$tip_id), , drop = FALSE]
}
