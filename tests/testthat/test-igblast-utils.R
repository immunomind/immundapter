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
