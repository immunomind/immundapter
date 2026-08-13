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
      make_igblast_test_immundata(sequence = NA_character_),
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
