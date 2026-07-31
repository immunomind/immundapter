test_that("IgBLAST runs ImmunData through preserved AIRR artifacts", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = function() {
      list(
        package_installed = TRUE,
        igblastr_version = "1.0.23",
        igblast_installed = TRUE,
        igblast_version = "1.22.0",
        igblast_root = "/mock/igblast",
        germline_db = "IMGT.mock.Homo_sapiens.TRB",
        c_region_db = ""
      )
    },
    igblast_call = function(arguments) {
      fasta <- readLines(arguments$query)
      query_ids <- sub("^>", "", fasta[seq.int(1L, length(fasta), by = 2L)])
      sequences <- fasta[seq.int(2L, length(fasta), by = 2L)]

      utils::write.table(
        data.frame(
          sequence_id = query_ids,
          sequence = sequences,
          locus = "TRB",
          v_call = "TRBV10-3*01",
          d_call = "TRBD1*01",
          j_call = "TRBJ2-3*01",
          productive = TRUE,
          sequence_alignment = sequences,
          germline_alignment = gsub("N", "A", sequences)
        ),
        file = arguments$out,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        na = ""
      )
      invisible(NULL)
    },
    .package = "immundapter"
  )

  idata <- make_igblast_test_immundata()
  run <- tool_igblast(conf_igblast(
    num_threads = 3L,
    run_args = list(num_alignments_V = 1L)
  ))$run(
    idata,
    output_dir = withr::local_tempdir()
  )

  expect_s3_class(run, "IgblastRun")
  expect_equal(run$tool, "igblast")
  expect_equal(run$config$num_threads, 3L)
  expect_setequal(
    names(run$artifacts),
    c("input", "config", "mapping", "airr")
  )
  expect_true(all(vapply(
    run$artifacts,
    fs::file_exists,
    logical(1L)
  )))
  expect_equal(
    readLines(run$artifacts$input),
    c(">igblast_query_000000001", "ACGTRYSWKMBDHVN")
  )

  manifest <- jsonlite::read_json(
    run$artifacts$config,
    simplifyVector = TRUE
  )
  expect_equal(manifest$fixed$outfmt, "AIRR")
  expect_equal(manifest$fixed$organism, "auto")
  expect_equal(manifest$fixed$ig_seqtype, "auto")
  expect_equal(manifest$fixed$domain_system, "imgt")
  expect_equal(
    manifest$runtime$germline_db,
    "IMGT.mock.Homo_sapiens.TRB"
  )

  results <- run$results()
  expect_s3_class(results, "IgblastResults")
  expect_identical(results$run, run)
  expect_s3_class(results$airr, "duckplyr_df")
  expect_s3_class(results$mapping, "duckplyr_df")
  expect_identical(results$raw_output, run$artifacts$airr)

  airr <- results$airr |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(airr), 1L)
  expect_equal(airr$sequence_id, "igblast_query_000000001")
  expect_equal(airr$v_call, "TRBV10-3*01")

  mapping <- results$mapping |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(mapping), 2L)
  expect_setequal(mapping$imd_chain_id, c("chain-1", "chain-2"))

  annotated <- results$annotate(idata)
  annotated_alignments <- annotated$annotations |>
    dplyr::select(
      imd_chain_id,
      sequence_alignment,
      germline_alignment
    ) |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(annotated_alignments), 3L)
  expect_true(all(!is.na(
    annotated_alignments$sequence_alignment[
      annotated_alignments$imd_chain_id %in% c("chain-1", "chain-2")
    ]
  )))
  expect_true(is.na(
    annotated_alignments$sequence_alignment[
      annotated_alignments$imd_chain_id == "chain-3"
    ]
  ))

  custom <- results$annotate(
    idata,
    columns = c("v_call", "productive")
  )
  expect_true(all(
    c("v_call", "productive") %in% colnames(custom$annotations)
  ))
  expect_error(
    results$annotate(idata, columns = "not_an_airr_column"),
    "result column(s) not found: not_an_airr_column",
    fixed = TRUE
  )
  expect_error(
    results$annotate(idata, columns = "sequence_id"),
    "cannot be annotated",
    fixed = TRUE
  )
})

test_that("IgBLAST accepts only ImmunData input", {
  input <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">query", "ACGT"), input)

  expect_error(
    tool_igblast()$run(
      input,
      output_dir = withr::local_tempdir()
    ),
    "Must be an R6 class"
  )
})

test_that("IgBLAST runs the human TCR fixture with igblastr", {
  testthat::skip_on_cran()
  testthat::skip_if_not_installed("igblastr")
  testthat::skip_if_not(igblastr::has_igblast())

  databases <- igblastr::list_germline_dbs(names.only = TRUE)
  database <- grep(
    "IMGT-.*Homo_sapiens.*TRA\\+TRB\\+TRG\\+TRD$",
    databases,
    value = TRUE
  )
  testthat::skip_if(
    length(database) == 0L,
    "A human TCR germline database is not installed."
  )
  suppressMessages(igblastr::use_germline_db(database[[1L]]))

  run <- tool_igblast()$run(
    read_igblast_fixture(),
    output_dir = withr::local_tempdir()
  )
  results <- run$results()
  airr <- results$airr |>
    dplyr::select(sequence_id, locus, v_call, j_call) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(airr), 3L)
  expect_true(all(airr$locus == "TRB"))
  expect_true(all(grepl("^TRBV", airr$v_call)))
  expect_true(all(grepl("^TRBJ", airr$j_call)))

  annotated <- results$annotate(read_igblast_fixture())
  expect_true(all(
    c("sequence_alignment", "germline_alignment") %in%
      colnames(annotated$annotations)
  ))
})
