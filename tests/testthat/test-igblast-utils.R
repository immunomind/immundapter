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
    "absent from the IgBLAST input: unknown",
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
  final <- igblast_merge_annotations(
    airr,
    annotations,
    "contig_id"
  ) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(final), 2L)
  expect_true(is.na(final$v_call[final$sequence_id == "query-2"]))
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

  prepared <- igblast_prepare_input(input)

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

  prepared <- igblast_prepare_input(input)

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
    c(fasta = fasta, annotations = annotations)
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
      c(fasta = fasta, annotations = annotations)
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
      c(fasta = fasta, annotations = annotations)
    ),
    "must be unique; duplicated: contig-1",
    fixed = TRUE
  )
})

test_that("IgBLAST validates FASTA structure and nucleotide sequences", {
  missing_header <- withr::local_tempfile(fileext = ".fasta")
  writeLines("ACGT", missing_header)
  expect_error(
    igblast_prepare_input(missing_header),
    "must start with a `>` header",
    fixed = TRUE
  )

  duplicate_id <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query", "ACGT", ">query description", "TGCA"), duplicate_id)
  expect_error(
    igblast_prepare_input(duplicate_id),
    "sequence names must be unique",
    fixed = TRUE
  )

  missing_sequence <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query-1", ">query-2", "ACGT"), missing_sequence)
  expect_error(
    igblast_prepare_input(missing_sequence),
    "contain no sequence: query-1",
    fixed = TRUE
  )

  invalid_sequence <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query", "CASSLG"), invalid_sequence)
  expect_error(
    igblast_prepare_input(invalid_sequence),
    "ungapped IUPAC nucleotide sequences",
    fixed = TRUE
  )
})

test_that("IgBLAST fixture contains full human TRB nucleotide sequences", {
  fixture <- read_igblast_fixture() |>
    dplyr::select(source_id, locus, sequence) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(fixture), 3L)
  expect_true(all(fixture$locus == "TRB"))
  expect_true(all(nchar(fixture$sequence) == 248L))
})

test_that("IMGT V lookup maps ungapped end-coordinate boundaries", {
  plan <- igblast_make_imgt_v_plan("synthetic*01", "AA..CC...GG")

  expect_equal(
    plan$imgt_v_end[match(2:5, plan$raw_v_end)],
    c(2L, 5L, 6L, 10L)
  )
  expect_equal(
    is.na(plan$pattern[match(2:5, plan$raw_v_end)]),
    c(TRUE, FALSE, FALSE, FALSE)
  )
})
