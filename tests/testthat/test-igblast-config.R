test_that("IgBLAST exposes only adapter-owned user settings", {
  expect_identical(
    names(formals(conf_igblast)),
    c("sequence_col", "num_threads", "run_args")
  )

  config <- conf_igblast()
  expect_equal(config$sequence_col, "sequence")
  expect_null(config$num_threads)
  expect_equal(config$run_args, list())
  expect_false("organism" %in% names(config))
  expect_false("ig_seqtype" %in% names(config))
  expect_false("igblast_root" %in% names(config))

  threaded <- conf_igblast(
    num_threads = 4L,
    run_args = list(num_alignments_V = 1L)
  )
  expect_equal(threaded$num_threads, 4L)
  expect_equal(threaded$run_args$num_alignments_V, 1L)
})

test_that("IgBLAST run arguments cannot override fixed behavior", {
  expect_error(
    conf_igblast(run_args = list(organism = "human")),
    "cannot override: organism",
    fixed = TRUE
  )
  expect_error(
    conf_igblast(run_args = list(domain_system = "kabat")),
    "cannot override: domain_system",
    fixed = TRUE
  )
  expect_error(
    conf_igblast(run_args = list(num_threads = 8L)),
    "cannot override: num_threads",
    fixed = TRUE
  )
  expect_error(
    conf_igblast(run_args = list(1L)),
    "unique, non-empty names",
    fixed = TRUE
  )
})

test_that("IgBLAST help uses igblastr and pak without installing", {
  output <- paste(
    testthat::capture_messages(tool_igblast()$help()),
    collapse = "\n"
  )

  expect_match(output, 'pak::pak\\("bioc::igblastr"\\)')
  expect_match(output, "igblastr::install_igblast")
  expect_match(output, "igblastr::install_IMGT_germline_db")
  expect_match(output, "igblastr::use_germline_db")
  expect_match(output, "C-region database is optional")
})

test_that("IgBLAST check reports a missing package without installing it", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = function() {
      list(
        package_installed = FALSE,
        igblast_installed = FALSE
      )
    },
    .package = "immundapter"
  )

  output <- testthat::capture_messages({
    checked <- tool_igblast()$check()
  })
  output <- paste(output, collapse = "\n")

  expect_false(checked)
  expect_match(output, "igblastr.*not installed")
})
