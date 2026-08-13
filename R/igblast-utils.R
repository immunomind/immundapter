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


igblast_prepare_immundata <- function(idata, sequence_col) {
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


igblast_read_fasta <- function(path) {
  connection <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, open = "rt")
  } else {
    file(path, open = "rt")
  }
  on.exit(close(connection), add = TRUE)
  lines <- readLines(connection, warn = FALSE)

  if (length(lines) > 0L) {
    lines[[1L]] <- sub("^\\ufeff", "", lines[[1L]])
  }
  lines <- lines[nzchar(trimws(lines)) & !grepl("^;", lines)]
  if (length(lines) == 0L) {
    rlang::abort("IgBLAST FASTA input contains no records.")
  }

  header_rows <- which(startsWith(lines, ">"))
  if (length(header_rows) == 0L || header_rows[[1L]] != 1L) {
    rlang::abort(
      "IgBLAST FASTA input must start with a `>` header line."
    )
  }

  headers <- trimws(sub("^>", "", lines[header_rows]))
  if (any(headers == "")) {
    rlang::abort("IgBLAST FASTA headers must not be empty.")
  }
  fasta_ids <- sub("[[:space:]].*$", "", headers)
  duplicate_ids <- unique(fasta_ids[duplicated(fasta_ids)])
  if (length(duplicate_ids) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST FASTA sequence names must be unique; repeated: %s.",
      paste(duplicate_ids, collapse = ", ")
    ))
  }

  end_rows <- c(header_rows[-1L] - 1L, length(lines))
  sequences <- vapply(seq_along(header_rows), function(index) {
    first_sequence_row <- header_rows[[index]] + 1L
    if (first_sequence_row > end_rows[[index]]) {
      return("")
    }
    sequence_rows <- seq.int(first_sequence_row, end_rows[[index]])
    gsub(
      "[[:space:]]",
      "",
      paste0(lines[sequence_rows], collapse = "")
    )
  }, character(1L))

  empty <- fasta_ids[sequences == ""]
  if (length(empty) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST FASTA record(s) contain no sequence: %s.",
      paste(empty, collapse = ", ")
    ))
  }

  data.frame(
    fasta_id = fasta_ids,
    fasta_header = headers,
    sequence = sequences,
    stringsAsFactors = FALSE
  )
}


igblast_prepare_fasta <- function(path) {
  records <- igblast_read_fasta(path) |>
    duckplyr::as_duckdb_tibble() |>
    dplyr::mutate(sequence = toupper(sequence))

  queries <- records |>
    dplyr::transmute(
      igblast_query_id = !!rlang::sym("fasta_id"),
      sequence = !!rlang::sym("sequence")
    ) |>
    dplyr::collect() |>
    as.data.frame()
  if (any(!grepl("^[ACGTRYSWKMBDHVN]+$", queries$sequence))) {
    rlang::abort(
      "IgBLAST FASTA records must contain ungapped IUPAC nucleotide sequences."
    )
  }

  mapping <- records |>
    dplyr::mutate(igblast_query_id = !!rlang::sym("fasta_id")) |>
    dplyr::select(dplyr::all_of(c(
      "fasta_id",
      "fasta_header",
      "igblast_query_id"
    ))) |>
    dplyr::distinct()

  list(queries = queries, mapping = mapping)
}


igblast_prepare_fasta_annotations <- function(fasta_path, annotations_path) {
  prepared <- igblast_prepare_fasta(fasta_path)
  delimiter <- if (
    grepl("\\.(tsv|txt)(\\.gz)?$", annotations_path, ignore.case = TRUE)
  ) {
    "\t"
  } else if (
    grepl("\\.csv(\\.gz)?$", annotations_path, ignore.case = TRUE)
  ) {
    ","
  } else {
    rlang::abort(
      "IgBLAST annotations must be a CSV or TSV file, optionally gzip-compressed."
    )
  }
  annotations <- duckplyr::read_csv_duckdb(
    annotations_path,
    options = list(delim = delimiter, header = TRUE)
  )
  annotation_id_col <- if ("contig_id" %in% colnames(annotations)) {
    "contig_id"
  } else if ("sequence_id" %in% colnames(annotations)) {
    "sequence_id"
  } else {
    rlang::abort(
      "IgBLAST annotations must contain `contig_id` or `sequence_id`."
    )
  }

  annotation_ids <- annotations |>
    dplyr::transmute(
      source_id = !!rlang::sym(annotation_id_col)
    )
  invalid_ids <- annotation_ids |>
    dplyr::filter(
      is.na(!!rlang::sym("source_id")) |
        !!rlang::sym("source_id") == ""
    ) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(invalid_ids) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST annotation `%s` values must not be missing or empty.",
      annotation_id_col
    ))
  }

  duplicate_ids <- annotation_ids |>
    dplyr::count(!!rlang::sym("source_id"), name = "id_count") |>
    dplyr::filter(!!rlang::sym("id_count") > 1L) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(duplicate_ids) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST annotation `%s` values must be unique; duplicated: %s.",
      annotation_id_col,
      paste(utils::head(duplicate_ids$source_id, 10L), collapse = ", ")
    ))
  }

  fasta_ids <- duckplyr::as_duckdb_tibble(data.frame(
    source_id = prepared$queries$igblast_query_id,
    stringsAsFactors = FALSE
  ))
  fasta_only <- fasta_ids |>
    dplyr::anti_join(annotation_ids, by = "source_id") |>
    dplyr::collect() |>
    as.data.frame()
  annotations_only <- annotation_ids |>
    dplyr::anti_join(fasta_ids, by = "source_id") |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(fasta_only) > 0L || nrow(annotations_only) > 0L) {
    rlang::abort(sprintf(
      paste(
        "IgBLAST FASTA sequence names and annotation-table IDs must match;",
        "%d occur only in FASTA and %d occur only in the annotation table."
      ),
      nrow(fasta_only),
      nrow(annotations_only)
    ))
  }

  prepared$annotations <- annotations
  prepared$annotation_id_col <- annotation_id_col
  prepared
}


igblast_merge_annotations <- function(airr, annotations, annotation_id_col) {
  airr_columns <- colnames(airr)
  annotation_columns <- colnames(annotations)
  if (!"sequence_id" %in% airr_columns) {
    rlang::abort("IgBLAST AIRR output is missing `sequence_id`.")
  }

  semantic_replacements <- c(
    chain = "locus",
    v_gene = "v_call",
    d_gene = "d_call",
    j_gene = "j_call",
    c_gene = "c_call",
    fwr1_nt = "fwr1",
    cdr1_nt = "cdr1",
    fwr2_nt = "fwr2",
    cdr2_nt = "cdr2",
    fwr3_nt = "fwr3",
    cdr3_nt = "junction",
    fwr4_nt = "fwr4"
  )
  semantic_collisions <- names(semantic_replacements)[
    names(semantic_replacements) %in% annotation_columns &
      unname(semantic_replacements) %in% airr_columns
  ]
  exact_collisions <- intersect(
    setdiff(annotation_columns, annotation_id_col),
    airr_columns
  )
  old_columns <- unique(c(exact_collisions, semantic_collisions))
  old_names <- paste0("old_", old_columns)
  occupied <- union(setdiff(annotation_columns, old_columns), airr_columns)
  conflicting_old_names <- intersect(old_names, occupied)
  if (length(conflicting_old_names) > 0L) {
    rlang::abort(sprintf(
      "Cannot preserve source annotations because column(s) already exist: %s.",
      paste(conflicting_old_names, collapse = ", ")
    ))
  }

  if (length(old_columns) > 0L) {
    rename_map <- stats::setNames(old_columns, old_names)
    annotations <- annotations |>
      dplyr::rename(dplyr::all_of(rename_map))
  }

  if (annotation_id_col == "sequence_id") {
    merged <- annotations |>
      dplyr::left_join(airr, by = "sequence_id")
  } else {
    merged <- annotations |>
      dplyr::mutate(
        sequence_id = !!rlang::sym(annotation_id_col)
      ) |>
      dplyr::left_join(airr, by = "sequence_id")
  }

  merged |>
    dplyr::relocate(dplyr::all_of("sequence_id"), .before = 1L)
}


igblast_prepare_input <- function(input, sequence_col) {
  if (inherits(input, "ImmunData")) {
    return(igblast_prepare_immundata(input, sequence_col))
  }
  if (length(input) == 2L) {
    return(igblast_prepare_fasta_annotations(
      input[["fasta"]],
      input[["annotations"]]
    ))
  }
  igblast_prepare_fasta(input)
}


igblast_call <- function(arguments) {
  invisible(do.call(
    getExportedValue("igblastr", "igblastn"),
    arguments
  ))
}
