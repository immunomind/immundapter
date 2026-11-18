test_that("RAxML-NG accepts only ImmunData input", {
  input <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">tip", "ACGT"), input)

  expect_error(
    tool_raxml_ng()$run(input, output_dir = withr::local_tempdir()),
    "ImmunData",
    fixed = TRUE
  )
})


test_that("RAxML-NG prepares ImmunData once per run", {
  unquote_arg <- function(value) {
    sub("^['\"](.*)['\"]$", "\\1", value)
  }
  original_prepare <- raxml_ng_prepare_immundata
  prepare_calls <- 0L
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_prepare_immundata = function(input, config) {
      prepare_calls <<- prepare_calls + 1L
      original_prepare(input, config)
    },
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      prefix <- unquote_arg(args[match("--prefix", args) + 1L])
      writeLines("(tip_000001,tip_000002);", paste0(
        prefix,
        ".raxml.bestTree"
      ))
      writeLines(character(), stdout)
      writeLines(character(), stderr)
      0L
    },
    .package = "immundapter"
  )

  tool_raxml_ng()$run(
    make_raxml_test_idata(),
    output_dir = withr::local_tempdir()
  )

  expect_equal(prepare_calls, 1L)
})


test_that("RAxML-NG includes additional tips in backend artifacts", {
  unquote_arg <- function(value) {
    sub("^['\"](.*)['\"]$", "\\1", value)
  }
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      prefix <- unquote_arg(args[match("--prefix", args) + 1L])
      msa <- unquote_arg(args[match("--msa", args) + 1L])
      tips <- sub("^>", "", grep("^>", readLines(msa), value = TRUE))
      writeLines(
        paste0("(", paste(tips, collapse = ","), ");"),
        paste0(prefix, ".raxml.bestTree")
      )
      writeLines(character(), stdout)
      writeLines(character(), stderr)
      0L
    },
    .package = "immundapter"
  )

  idata <- make_raxml_test_idata()
  additional_tips <- raxml_test_additional_tips() |>
    dplyr::filter(family_id == "family-a")

  run <- tool_raxml_ng()$run(
    idata,
    additional_tips = additional_tips,
    output_dir = withr::local_tempdir()
  )
  normalized <- utils::read.delim(
    run$artifacts$input,
    stringsAsFactors = FALSE
  )
  tip_map <- run$results()$tip_map |>
    dplyr::collect() |>
    as.data.frame()
  outgroups <- run$results()$outgroups |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(run$status, "success")
  expect_equal(nrow(normalized), 3L)
  expect_true("family-a__germline" %in% normalized$tip_id)
  expect_true("family-a__germline" %in% tip_map$tip_id)
  expect_equal(nrow(outgroups), 0L)
  expect_equal(
    idata$annotations |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::collect() |>
      dplyr::pull(n),
    2L
  )
})


test_that("RAxML-NG roots family trees and reconstructs ancestral states", {
  unquote_arg <- function(value) {
    sub("^['\"](.*)['\"]$", "\\1", value)
  }
  calls <- list()
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      calls[[length(calls) + 1L]] <<- args
      prefix <- unquote_arg(args[match("--prefix", args) + 1L])
      msa <- unquote_arg(args[match("--msa", args) + 1L])
      tips <- sub("^>", "", grep("^>", readLines(msa), value = TRUE))
      tree <- paste0(
        "(",
        paste0(tips, ":0.1", collapse = ","),
        ")Node1;"
      )
      writeLines(character(), stdout)
      writeLines(character(), stderr)

      if ("--ancestral" %in% args) {
        writeLines(tree, paste0(prefix, ".raxml.ancestralTree"))
        writeLines(
          "Node1\tACGT",
          paste0(prefix, ".raxml.ancestralStates")
        )
        writeLines(
          c(
            "Node\tSite\tState\tp_A\tp_C\tp_G\tp_T",
            "Node1\t1\tA\t1\t0\t0\t0"
          ),
          paste0(prefix, ".raxml.ancestralProbs")
        )
      } else {
        writeLines(tree, paste0(prefix, ".raxml.bestTree"))
      }
      writeLines("mock success", paste0(prefix, ".raxml.log"))
      0L
    },
    .package = "immundapter"
  )

  additional_tips <- raxml_test_additional_tips() |>
    dplyr::filter(family_id == "family-a")
  outgroups <- additional_tips |>
    dplyr::select(family_id, tip_id)
  run <- tool_raxml_ng(conf_raxml_ng(ancestral = TRUE))$run(
    make_raxml_test_idata(),
    additional_tips = additional_tips,
    outgroups = outgroups,
    output_dir = withr::local_tempdir()
  )

  results <- run$results()
  tree_index <- results$tree_index |>
    dplyr::collect() |>
    as.data.frame()
  ancestral_index <- results$ancestral_index |>
    dplyr::collect() |>
    as.data.frame()
  outgroup_map <- results$outgroups |>
    dplyr::collect() |>
    as.data.frame()
  states <- results$read_ancestral_states() |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(run$status, "success")
  expect_equal(length(calls), 2L)
  expect_true("--ancestral" %in% calls[[2L]])
  expect_false(any(vapply(calls, function(args) {
    "--outgroup" %in% args
  }, logical(1L))))
  expect_equal(outgroup_map$tip_id, "family-a__germline")
  expect_false(ape::is.rooted(ape::read.tree(tree_index$best_tree)))
  expect_true(ape::is.rooted(ape::read.tree(tree_index$tree)))
  expect_true(ape::is.rooted(ape::read.tree(ancestral_index$ancestral_tree)))
  expect_equal(states$family_id, "family-a")
  expect_equal(states$node_id, "Node1")
  expect_equal(states$sequence, "ACGT")
  expect_named(results$read_ancestral_trees(), "family-a")
})


test_that("RAxML-NG preserves a tree when ancestral reconstruction fails", {
  unquote_arg <- function(value) {
    sub("^['\"](.*)['\"]$", "\\1", value)
  }
  testthat::local_mocked_bindings(
    raxml_ng_executable = function() "/mock/raxml-ng",
    raxml_ng_execute = function(executable, args, stdout, stderr) {
      prefix <- unquote_arg(args[match("--prefix", args) + 1L])
      writeLines(character(), stdout)
      if ("--ancestral" %in% args) {
        writeLines("mock ancestral failure", stderr)
        return(8L)
      }
      writeLines("(tip_000001:0.1,tip_000002:0.1);", paste0(
        prefix,
        ".raxml.bestTree"
      ))
      writeLines(character(), stderr)
      0L
    },
    .package = "immundapter"
  )

  run <- tool_raxml_ng(conf_raxml_ng(ancestral = TRUE))$run(
    make_raxml_test_idata(),
    output_dir = withr::local_tempdir()
  )
  results <- run$results()
  failed <- results$failed_ancestral |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(run$status, "partial")
  expect_named(results$read_trees(), "family-a")
  expect_equal(failed$family_id, "family-a")
  expect_match(failed$ancestral_error, "mock ancestral failure", fixed = TRUE)
  expect_length(results$read_ancestral_trees(), 0L)
  expect_equal(
    results$read_ancestral_states() |>
      dplyr::summarise(n = dplyr::n()) |>
      dplyr::collect() |>
      dplyr::pull(n),
    0L
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
  ancestral_run <- tool_raxml_ng(conf_raxml_ng(
    analysis = "search",
    ancestral = TRUE,
    seed = 23,
    threads = 1
  ))$run(
    idata,
    outgroups = data.frame(
      family_id = "family-a",
      tip_id = "1"
    ) |>
      duckplyr::as_duckdb_tibble(),
    output_dir = output_dir
  )

  expect_true(inherits(run, "RaxmlNgRun"))
  expect_true(inherits(fast_run, "RaxmlNgRun"))
  expect_true(inherits(ancestral_run, "RaxmlNgRun"))
  expect_equal(run$status, "success")
  expect_equal(ancestral_run$status, "success")
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

  ancestral_results <- ancestral_run$results()
  ancestral_index <- ancestral_results$ancestral_index |>
    dplyr::collect() |>
    as.data.frame()
  rooted_index <- ancestral_results$tree_index |>
    dplyr::collect() |>
    as.data.frame()
  ancestral_states <- ancestral_results$read_ancestral_states() |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(ancestral_index), 1L)
  expect_true(ape::is.rooted(ape::read.tree(rooted_index$tree)))
  expect_true(ape::is.rooted(ape::read.tree(
    ancestral_index$ancestral_tree
  )))
  expect_gt(nrow(ancestral_states), 0L)
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

  idata <- make_raxml_test_idata()
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
