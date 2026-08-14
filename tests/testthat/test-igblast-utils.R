test_that("IgBLAST prepares and deduplicates ImmunData chain sequences", {
  idata <- make_igblast_test_immundata()
  prepared <- igblast_prepare_input(idata, "sequence")

  expect_equal(
    prepared$queries,
    data.frame(
      igblast_query_id = "igblast_query_000000001",
      sequence = "ACGTRYSWKMBDHVN"
    )
  )
  expect_s3_class(prepared$mapping, "duckplyr_df")

  mapping <- prepared$mapping |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(mapping), 2L)
  expect_setequal(mapping$imd_chain_id, c("chain-1", "chain-2"))
  expect_equal(length(unique(mapping$igblast_query_id)), 1L)
})

test_that("IgBLAST retains missing and empty ImmunData input status", {
  prepared <- igblast_prepare_input(
    make_igblast_test_immundata(
      sequence = c("ACGT", NA_character_, ""),
      chain_id = c("chain-submitted", "chain-missing", "chain-empty")
    ),
    "sequence"
  )
  input_status <- prepared$input_status |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(input_status), 3L)
  expect_equal(
    input_status$submission_status[
      input_status$input_id == "chain-submitted"
    ],
    "submitted"
  )
  expect_equal(
    input_status$submission_status[input_status$input_id == "chain-missing"],
    "excluded_missing_sequence"
  )
  expect_equal(
    input_status$submission_status[input_status$input_id == "chain-empty"],
    "excluded_empty_sequence"
  )
  expect_true(is.na(
    input_status$igblast_query_id[input_status$input_id == "chain-missing"]
  ))
  expect_true(is.na(
    input_status$igblast_query_id[input_status$input_id == "chain-empty"]
  ))
})

test_that("IgBLAST validates AIRR IDs and classifies every input", {
  submitted <- igblast_test_submitted_queries()
  airr <- duckplyr::as_duckdb_tibble(data.frame(
    sequence_id = c("query-1", "query-2"),
    v_call = c("TRBV1*01", "TRBV2*01"),
    j_call = c("TRBJ1*01", NA_character_)
  ))
  assignments <- igblast_validate_airr(airr, submitted)
  input_status <- duckplyr::as_duckdb_tibble(data.frame(
    input_id = c("chain-1", "chain-2", "chain-3", "chain-4", "chain-5"),
    input_id_type = "imd_chain_id",
    igblast_query_id = c(
      "query-1", "query-1", "query-2", "query-3", NA_character_
    ),
    submission_status = c(
      "submitted", "submitted", "submitted", "submitted",
      "excluded_empty_sequence"
    ),
    submission_reason = c(
      NA_character_, NA_character_, NA_character_, NA_character_,
      "empty_sequence"
    )
  ))

  status <- igblast_make_query_status(input_status, assignments) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(
    status$result_status,
    c("assigned", "assigned", "unassigned", "missing_output", "not_submitted")
  )
  expect_equal(
    status$reason,
    c(
      NA_character_, NA_character_, "missing_j_call",
      "absent_from_airr_output", "empty_sequence"
    )
  )
})

test_that("IgBLAST rejects unknown, duplicate, and malformed AIRR IDs", {
  submitted <- igblast_test_submitted_queries()[1:2, ]

  expect_error(
    igblast_validate_airr(
      duckplyr::as_duckdb_tibble(data.frame(
        sequence_id = "unknown",
        v_call = "TRBV1*01",
        j_call = "TRBJ1*01"
      )),
      submitted
    ),
    "were not submitted: unknown",
    fixed = TRUE
  )
  expect_error(
    igblast_validate_airr(
      duckplyr::as_duckdb_tibble(data.frame(
        sequence_id = c("query-1", "query-1"),
        v_call = "TRBV1*01",
        j_call = "TRBJ1*01"
      )),
      submitted
    ),
    "must be unique; duplicated: query-1",
    fixed = TRUE
  )
  expect_error(
    igblast_validate_airr(
      duckplyr::as_duckdb_tibble(data.frame(
        sequence_id = c("query-1", " "),
        v_call = "TRBV1*01",
        j_call = "TRBJ1*01"
      )),
      submitted
    ),
    "must not be missing or empty",
    fixed = TRUE
  )
  expect_error(
    igblast_validate_airr(
      duckplyr::as_duckdb_tibble(data.frame(v_call = "TRBV1*01")),
      submitted
    ),
    "missing `sequence_id`",
    fixed = TRUE
  )
})

test_that("IgBLAST paired merge preserves rows missing from AIRR output", {
  airr <- duckplyr::as_duckdb_tibble(data.frame(
    sequence_id = "query-1",
    v_call = "IGHV1*01",
    j_call = "IGHJ1*01"
  ))
  annotations <- duckplyr::as_duckdb_tibble(data.frame(
    contig_id = c("query-1", "query-2"),
    barcode = c("cell-1", "cell-2")
  ))
  query_status <- duckplyr::as_duckdb_tibble(data.frame(
    input_id = c("query-1", "query-2"),
    input_id_type = "fasta_id",
    igblast_query_id = c("query-1", "query-2"),
    submission_status = "submitted",
    result_status = c("assigned", "missing_output"),
    reason = c(NA_character_, "absent_from_airr_output")
  ))

  final <- igblast_merge_annotations(
    airr,
    annotations,
    "contig_id",
    query_status
  ) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(final), 2L)
  expect_equal(final$result_status, c("assigned", "missing_output"))
  expect_true(is.na(final$v_call[final$sequence_id == "query-2"]))
})

test_that("IgBLAST validates its ImmunData sequence input", {
  expect_error(
    igblast_prepare_input(
      make_igblast_test_immundata(sequence = "CASSLG"),
      "sequence"
    ),
    "ungapped IUPAC nucleotide",
    fixed = TRUE
  )
  expect_error(
    igblast_prepare_input(
      make_igblast_test_immundata(sequence = "ACGT-ACGT"),
      "sequence"
    ),
    "ungapped IUPAC nucleotide",
    fixed = TRUE
  )
  expect_error(
    igblast_prepare_input(
      make_igblast_test_immundata(
        sequence = c(NA_character_, ""),
        chain_id = c("chain-missing", "chain-empty")
      ),
      "sequence"
    ),
    "No nucleotide sequences",
    fixed = TRUE
  )
})

test_that("IgBLAST prepares multiline gzip-compressed FASTA input", {
  input <- withr::local_tempfile(fileext = ".fasta.gz")
  connection <- gzfile(input, open = "wt")
  writeLines(c(
    ">sequence-1 description",
    "acgt",
    "rysw",
    ">sequence-2",
    "TTAA"
  ), connection)
  close(connection)

  prepared <- igblast_prepare_input(input, "unused")

  expect_equal(
    prepared$queries,
    data.frame(
      igblast_query_id = c("sequence-1", "sequence-2"),
      sequence = c("ACGTRYSW", "TTAA")
    )
  )
  expect_s3_class(prepared$mapping, "duckplyr_df")
  mapping <- prepared$mapping |>
    dplyr::collect() |>
    as.data.frame()
  expect_setequal(mapping$fasta_id, c("sequence-1", "sequence-2"))
  expect_equal(
    mapping$fasta_header[mapping$fasta_id == "sequence-1"],
    "sequence-1 description"
  )
  expect_equal(mapping$igblast_query_id, mapping$fasta_id)
})

test_that("IgBLAST preserves FASTA sequence names and duplicate sequences", {
  input <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(
    ">contig-2 description",
    "acgt",
    ">contig-1",
    "ACGT",
    ">contig-3",
    "ttgg"
  ), input)

  prepared <- igblast_prepare_input(input, "unused")

  expect_equal(
    prepared$queries,
    data.frame(
      igblast_query_id = c("contig-2", "contig-1", "contig-3"),
      sequence = c("ACGT", "ACGT", "TTGG")
    )
  )
  expect_equal(nrow(prepared$queries), 3L)
})

test_that("IgBLAST validates paired FASTA and annotation-table names", {
  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(
    ">contig-1",
    "ACGT",
    ">contig-2",
    "TGCA"
  ), fasta)
  annotations <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      contig_id = c("contig-1", "contig-2"),
      barcode = c("cell-1", "cell-2")
    ),
    annotations,
    row.names = FALSE
  )

  prepared <- igblast_prepare_input(
    c(fasta = fasta, annotations = annotations),
    "unused"
  )
  expect_equal(prepared$annotation_id_col, "contig_id")
  expect_s3_class(prepared$annotations, "duckplyr_df")

  utils::write.csv(
    data.frame(contig_id = c("contig-1", "contig-3")),
    annotations,
    row.names = FALSE
  )
  expect_error(
    igblast_prepare_input(
      c(fasta = fasta, annotations = annotations),
      "unused"
    ),
    "1 occur only in FASTA and 1 occur only in the annotation table",
    fixed = TRUE
  )

  utils::write.csv(
    data.frame(contig_id = c("contig-1", "contig-1")),
    annotations,
    row.names = FALSE
  )
  expect_error(
    igblast_prepare_input(
      c(fasta = fasta, annotations = annotations),
      "unused"
    ),
    "must be unique; duplicated: contig-1",
    fixed = TRUE
  )
})

test_that("IgBLAST validates FASTA structure and nucleotide sequences", {
  missing_header <- withr::local_tempfile(fileext = ".fasta")
  writeLines("ACGT", missing_header)
  expect_error(
    igblast_prepare_input(missing_header, "unused"),
    "must start with a `>` header",
    fixed = TRUE
  )

  duplicate_id <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query", "ACGT", ">query description", "TGCA"), duplicate_id)
  expect_error(
    igblast_prepare_input(duplicate_id, "unused"),
    "sequence names must be unique",
    fixed = TRUE
  )

  missing_sequence <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query-1", ">query-2", "ACGT"), missing_sequence)
  expect_error(
    igblast_prepare_input(missing_sequence, "unused"),
    "contain no sequence: query-1",
    fixed = TRUE
  )

  invalid_sequence <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query", "CASSLG"), invalid_sequence)
  expect_error(
    igblast_prepare_input(invalid_sequence, "unused"),
    "ungapped IUPAC nucleotide sequences",
    fixed = TRUE
  )
})

test_that("IgBLAST fixture contains full human TRB nucleotide sequences", {
  idata <- read_igblast_fixture()
  fixture <- idata$annotations |>
    dplyr::select(source_id, locus, sequence) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(fixture), 3L)
  expect_true(all(fixture$locus == "TRB"))
  expect_true(all(nchar(fixture$sequence) == 248L))
})
