test_that("RAxML-NG accepts only ImmunData input", {
  input <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">tip", "ACGT"), input)

  expect_error(
    tool_raxml_ng()$run(input, output_dir = withr::local_tempdir()),
    "ImmunData",
    fixed = TRUE
  )
})


test_that("RAxML-NG infers and indexes family trees", {
  testthat::skip_on_cran()
  testthat::skip_if(Sys.which("raxml-ng") == "")

  reference <- paste(rep("ACGT", 30L), collapse = "")
  mutate_base <- function(sequence, position, base) {
    substring(sequence, position, position) <- base
    sequence
  }
  sequences <- c(
    reference,
    mutate_base(reference, 5L, "T"),
    mutate_base(reference, 17L, "C"),
    mutate_base(reference, 38L, "A"),
    mutate_base(reference, 71L, "A"),
    mutate_base(reference, 89L, "G")
  )
  annotations <- data.frame(
    imd_receptor_id = seq_along(sequences),
    clone_id = "family-a",
    clone_sequence_alignment = sequences,
    sequence = sequences
  ) |>
    duckplyr::as_duckdb_tibble()
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "sequence",
      chains = "IGH"
    ),
    annotations = annotations
  )

  output_dir <- withr::local_tempdir()
  run <- tool_raxml_ng(conf_raxml_ng(
    analysis = "search"
  ))$run(
    idata,
    output_dir = output_dir
  )
  fast_run <- tool_raxml_ng(conf_raxml_ng(
    analysis = "fast"
  ))$run(
    idata,
    output_dir = output_dir
  )

  expect_true(inherits(run, "RaxmlNgRun"))
  expect_true(inherits(fast_run, "RaxmlNgRun"))
  expect_equal(run$status, "success")
  expect_true(all(vapply(run$artifacts, fs::file_exists, logical(1L))))

  results <- run$results()
  expect_true(inherits(results, "RaxmlNgResults"))
  tree_index <- results$tree_index |>
    dplyr::collect() |>
    as.data.frame()
  tip_map <- results$tip_map |>
    dplyr::collect() |>
    as.data.frame()
  trees <- results$read_trees("family-a")

  expect_equal(nrow(tree_index), 1L)
  expect_equal(nrow(tip_map), length(sequences))
  expect_named(trees, "family-a")
  expect_match(unname(trees), "tip_000001", fixed = TRUE)
  expect_error(
    results$read_trees("missing-family"),
    "No successful RAxML-NG tree",
    fixed = TRUE
  )
  expect_error(
    results$annotate(idata),
    "cannot be added as columns to `ImmunData`",
    class = "immundapter_unannotatable_results"
  )
})
