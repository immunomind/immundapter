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
  runtime$igblast_build <- as.character(info$igblast_build)
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


igblast_find_bundled_c_region_db <- function(germline_db) {
  germline_dbs <- tryCatch(
    getExportedValue("igblastr", "list_germline_dbs")(
      long.listing = TRUE
    ),
    error = function(error) NULL
  )
  germline_loci <- rownames(germline_dbs[[germline_db]])
  if (is.null(germline_loci)) {
    return("")
  }

  organism_aliases <- list(
    human = c("human", "Homo_sapiens"),
    mouse = c("mouse", "Mus_musculus"),
    rabbit = c("rabbit", "Oryctolagus_cuniculus"),
    rat = c("rat", "Rattus_norvegicus"),
    rhesus_monkey = c("rhesus_monkey", "Macaca_mulatta")
  )
  organisms <- names(organism_aliases)[vapply(
    organism_aliases,
    function(aliases) any(vapply(
      aliases,
      grepl,
      logical(1),
      x = germline_db,
      fixed = TRUE
    )),
    logical(1)
  )]
  if (length(organisms) != 1L) {
    return("")
  }

  c_region_dbs <- tryCatch(
    getExportedValue("igblastr", "list_c_region_dbs")(
      builtin.only = TRUE,
      long.listing = TRUE
    ),
    error = function(error) NULL
  )
  candidates <- c_region_dbs[grepl(
    paste0(".", organisms, "."),
    names(c_region_dbs),
    fixed = TRUE
  )]
  if (length(candidates) == 0L) {
    return("")
  }

  overlap <- vapply(
    candidates,
    function(database) length(intersect(names(database), germline_loci)),
    integer(1)
  )
  if (max(overlap) == 0L || sum(overlap == max(overlap)) != 1L) {
    return("")
  }
  names(candidates)[which.max(overlap)]
}


igblast_make_imgt_v_plan <- function(v_allele, reference) {
  bases <- strsplit(reference, "", fixed = TRUE)[[1L]]
  ungapped_to_gapped <- which(bases != ".")
  runs <- rle(bases == ".")
  run_ends <- cumsum(runs$lengths)
  run_starts <- run_ends - runs$lengths + 1L
  gap_runs <- which(runs$values)
  gap_after <- vapply(
    run_starts[gap_runs],
    function(start) sum(bases[seq_len(start - 1L)] != "."),
    integer(1L)
  )
  gap_lengths <- runs$lengths[gap_runs]

  raw_v_end <- seq_along(ungapped_to_gapped)
  gap_count <- findInterval(raw_v_end - 1L, gap_after)
  patterns <- replacements <- rep(NA_character_, length(gap_lengths) + 1L)
  for (count in unique(gap_count[gap_count > 0L])) {
    chunks <- diff(c(0L, gap_after[seq_len(count)]))
    patterns[[count + 1L]] <- paste0(
      "^", paste0("(.{", chunks, "})", collapse = "")
    )
    replacements[[count + 1L]] <- paste0(
      paste0("\\", seq_len(count), strrep(".", gap_lengths[seq_len(count)])),
      collapse = ""
    )
  }

  data.frame(
    v_allele,
    raw_v_end = as.integer(raw_v_end),
    imgt_v_end = as.integer(ungapped_to_gapped),
    pattern = patterns[gap_count + 1L],
    replacement = replacements[gap_count + 1L]
  )
}


igblast_make_imgt_v_plans <- function(path) {
  records <- do.call(rbind, lapply(path, igblast_read_fasta))
  fields <- strsplit(records$fasta_header, "|", fixed = TRUE)
  v_alleles <- vapply(fields, function(parts) {
    if (length(parts) < 2L || !nzchar(parts[[2L]])) {
      rlang::abort("IMGT V FASTA headers must contain an allele in field 2.")
    }
    parts[[2L]]
  }, character(1L))
  sequence_counts <- vapply(
    split(records$sequence, v_alleles),
    function(sequence) length(unique(sequence)),
    integer(1L)
  )
  if (any(sequence_counts > 1L)) {
    rlang::abort("IMGT V FASTA alleles must not have conflicting sequences.")
  }
  keep <- !duplicated(v_alleles)
  records <- records[keep, , drop = FALSE]
  v_alleles <- v_alleles[keep]

  do.call(rbind, Map(igblast_make_imgt_v_plan, v_alleles, records$sequence))
}


igblast_find_imgt_v_references <- function(germline_db) {
  if (!startsWith(germline_db, "IMGT-")) {
    return(character())
  }

  directory <- file.path(
    tools::R_user_dir("igblastr", "cache"),
    "germline_dbs", germline_db, "V_original_fasta"
  )
  paths <- if (dir.exists(directory)) {
    list.files(directory, "[.]fasta$", full.names = TRUE)
  } else {
    character()
  }
  if (!length(paths)) {
    rlang::abort(paste(
      "The selected IMGT germline database does not contain its original",
      "gapped V FASTA files. Reinstall it with",
      "`igblastr::install_IMGT_germline_db()`."
    ))
  }
  paths
}


igblast_normalize_imgt <- function(airr, plans) {
  plan_table <- duckplyr::as_duckdb_tibble(plans, prudence = "stingy")
  input <- airr |>
    dplyr::mutate(
      .v_allele = dd$trim(dd$regexp_extract(
        v_call, "[^,|[:space:]][^,|]*", 0L
      ))
    )

  trimmed_insertions <- input |>
    dplyr::filter(dd$contains(germline_alignment, "-")) |>
    dplyr::mutate(
      .position = dd$unnest(dd$generate_series(
        1L, as.integer(dd$length(germline_alignment))
      )),
      .germline_base = dd$substr(germline_alignment, .position, 1L),
      .query_base = dplyr::if_else(
        .germline_base == "-",
        "",
        dd$substr(sequence_alignment, .position, 1L)
      )
    ) |>
    dplyr::arrange(sequence_id, .position) |>
    dplyr::summarise(
      .trimmed_sequence = dd$string_agg(.query_base, ""),
      .by = sequence_id
    )

  input |>
    dplyr::left_join(trimmed_insertions, by = "sequence_id") |>
    dplyr::mutate(
      .trimmed_sequence = dplyr::coalesce(
        .trimmed_sequence, sequence_alignment
      )
    ) |>
    dplyr::left_join(
      plan_table,
      by = c(".v_allele" = "v_allele", "v_germline_end" = "raw_v_end"),
      copy = TRUE
    ) |>
    dplyr::mutate(
      .normalized_sequence = dd$concat(
        dd$rpad("", as.integer(v_germline_start - 1L), "."),
        .trimmed_sequence
      ),
      .normalized_sequence = dplyr::if_else(
        is.na(pattern),
        .normalized_sequence,
        dd$regexp_replace(.normalized_sequence, pattern, replacement)
      ),
      sequence_alignment = dplyr::if_else(
        is.na(imgt_v_end), sequence_alignment, .normalized_sequence
      ),
      v_germline_start = dplyr::if_else(
        is.na(imgt_v_end), as.integer(v_germline_start), 1L
      ),
      v_germline_end = dplyr::coalesce(
        imgt_v_end, as.integer(v_germline_end)
      ),
      j_germline_start = as.integer(j_germline_start),
      j_germline_end = as.integer(j_germline_end)
    ) |>
    dplyr::select(-dplyr::all_of(c(
      ".v_allele", ".trimmed_sequence", ".normalized_sequence", "imgt_v_end",
      "pattern", "replacement"
    )))
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
    dplyr::mutate(sequence = dd$upper(sequence))

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
    )))

  list(
    queries = queries,
    mapping = mapping
  )
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
  annotation_header <- duckplyr::read_csv_duckdb(
    annotations_path,
    options = list(
      delim = delimiter,
      header = TRUE,
      all_varchar = TRUE
    )
  )
  annotation_id_col <- if ("contig_id" %in% colnames(annotation_header)) {
    "contig_id"
  } else if ("sequence_id" %in% colnames(annotation_header)) {
    "sequence_id"
  } else {
    rlang::abort(
      "IgBLAST annotations must contain `contig_id` or `sequence_id`."
    )
  }
  annotations <- duckplyr::read_csv_duckdb(
    annotations_path,
    options = list(
      delim = delimiter,
      header = TRUE,
      types = list(stats::setNames("VARCHAR", annotation_id_col))
    )
  )

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


igblast_merge_annotations <- function(
  airr,
  annotations,
  annotation_id_col
) {
  airr_columns <- colnames(airr)
  annotation_columns <- colnames(annotations)

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


igblast_validate_airr <- function(airr, queries) {
  if (!"sequence_id" %in% colnames(airr)) {
    rlang::abort("IgBLAST AIRR output is missing `sequence_id`.")
  }

  result_ids <- airr |>
    dplyr::select(dplyr::all_of("sequence_id"))
  invalid_ids <- result_ids |>
    dplyr::filter(
      is.na(!!rlang::sym("sequence_id")) |
        trimws(!!rlang::sym("sequence_id")) == ""
    ) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(invalid_ids) > 0L) {
    rlang::abort(
      "IgBLAST AIRR `sequence_id` values must not be missing or empty."
    )
  }

  duplicate_ids <- result_ids |>
    dplyr::count(!!rlang::sym("sequence_id"), name = "result_count") |>
    dplyr::filter(!!rlang::sym("result_count") > 1L) |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(duplicate_ids) > 0L) {
    rlang::abort(sprintf(
      "IgBLAST AIRR `sequence_id` values must be unique; duplicated: %s.",
      paste(utils::head(duplicate_ids$sequence_id, 10L), collapse = ", ")
    ))
  }

  input_ids <- duckplyr::as_duckdb_tibble(data.frame(
    sequence_id = queries$igblast_query_id,
    stringsAsFactors = FALSE
  ))
  unknown_ids <- result_ids |>
    dplyr::anti_join(input_ids, by = "sequence_id") |>
    dplyr::collect() |>
    as.data.frame()
  if (nrow(unknown_ids) > 0L) {
    rlang::abort(sprintf(
      paste(
        "IgBLAST AIRR output contains `sequence_id` values absent from",
        "the IgBLAST input: %s."
      ),
      paste(utils::head(unknown_ids$sequence_id, 10L), collapse = ", ")
    ))
  }

  invisible(NULL)
}


igblast_prepare_input <- function(input) {
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
