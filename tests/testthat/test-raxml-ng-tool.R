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


test_that("RAxML-NG returns successful families when another family fails", {
  unquote_arg <- function(value) {
    sub("^['\"](.*)['\"]$", "\\1", value)
  }
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      prefix <- unquote_arg(args[match("--prefix", args) + 1L])
      msa <- unquote_arg(args[match("--msa", args) + 1L])
      writeLines(character(), stdout)
      if (grepl("family-000002", prefix, fixed = TRUE)) {
        writeLines("mock family failure", stderr)
        return(7L)
      }

      tips <- sub("^>", "", grep("^>", readLines(msa), value = TRUE))
      writeLines(paste0("(", paste(tips, collapse = ","), ");"), paste0(
        prefix,
        ".raxml.bestTree"
      ))
      writeLines("mock success", paste0(prefix, ".raxml.log"))
      writeLines(character(), stderr)
      0L
    },
    .package = "immundapter"
  )

  annotations <- data.frame(
    imd_receptor_id = c("germline", "tip-a", "tip-b", "tip-c"),
    clone_id = rep(c("family-a", "family-b"), each = 2L),
    clone_sequence_alignment = c("ACGT", "ACGA", "TCGT", "TCGA"),
    partitions = "V=1-2;J=3-4",
    sequence = c("ACGT", "ACGA", "TCGT", "TCGA")
  ) |>
    duckplyr::as_duckdb_tibble()
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "sequence",
      chains = "IGH"
    ),
    annotations = annotations
  )
  run <- tool_raxml_ng(conf_raxml_ng(
    partition_col = "partitions",
    seed = 19,
    threads = 2
  ))$run(idata, output_dir = withr::local_tempdir())

  expect_s3_class(run, "RaxmlNgRun")
  expect_equal(run$status, "partial")
  expect_true("failed_families" %in% names(run$artifacts))
  expect_true(all(vapply(run$artifacts, fs::file_exists, logical(1L))))

  manifest <- jsonlite::read_json(
    run$artifacts$config,
    simplifyVector = TRUE
  )
  expect_equal(manifest$model, "GTR+G")
  expect_equal(manifest$seed, 19)
  expect_equal(manifest$threads, 2)
  expect_equal(manifest$brlen, "scaled")

  results <- run$results()
  families <- results$families |>
    dplyr::collect() |>
    as.data.frame()
  tree_index <- results$tree_index |>
    dplyr::collect() |>
    as.data.frame()
  failed <- results$failed_families |>
    dplyr::collect() |>
    as.data.frame()
  tip_map <- results$tip_map |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(families$status, c("success", "failed"))
  expect_true(all(fs::file_exists(families$model_file)))
  expect_true(all(grepl("--brlen scaled", families$command, fixed = TRUE)))
  expect_true(all(grepl("--seed 19", families$command, fixed = TRUE)))
  expect_true(all(grepl("--threads 2", families$command, fixed = TRUE)))
  expect_equal(families$exit_status, c(0L, 7L))
  expect_equal(tree_index$family_id, "family-a")
  expect_equal(failed$family_id, "family-b")
  expect_match(failed$error, "mock family failure", fixed = TRUE)
  expect_true("germline" %in% tip_map$tip_id)
  expect_named(results$read_trees(), "family-a")
  expect_error(
    results$read_trees("family-b"),
    "No successful RAxML-NG tree",
    fixed = TRUE
  )
})


test_that("RAxML-NG returns diagnostics when every family fails", {
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      writeLines(character(), stdout)
      writeLines("mock total failure", stderr)
      9L
    },
    .package = "immundapter"
  )

  annotations <- data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    clone_sequence_alignment = c("ACGT", "ACGA"),
    sequence = c("ACGT", "ACGA")
  ) |>
    duckplyr::as_duckdb_tibble()
  idata <- immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "sequence",
      chains = "IGH"
    ),
    annotations = annotations
  )
  run <- tool_raxml_ng()$run(
    idata,
    output_dir = withr::local_tempdir()
  )

  expect_equal(run$status, "partial")
  results <- run$results()
  families <- results$families |>
    dplyr::collect() |>
    as.data.frame()
  tree_index <- results$tree_index |>
    dplyr::collect() |>
    as.data.frame()
  failed <- results$failed_families |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(families$status, "failed")
  expect_equal(nrow(tree_index), 0L)
  expect_equal(nrow(failed), 1L)
  expect_match(failed$error, "mock total failure", fixed = TRUE)
  expect_length(results$read_trees(), 0L)
})
