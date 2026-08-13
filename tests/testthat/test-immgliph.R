make_immgliph_test_idata <- function() {
  annotations <- data.frame(
    imd_receptor_id = c("receptor-1", "receptor-2", "receptor-3"),
    locus = "TRB",
    junction_aa = c("CASSAAA", "CASSAAB", "CASSCCC"),
    v_call = c("TRBV1*01", "TRBV1*01", "TRBV2*01"),
    donor = c("donor-1", "donor-1", "donor-2"),
    stringsAsFactors = FALSE
  ) |>
    duckplyr::as_duckdb_tibble()

  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = c("junction_aa", "v_call"),
      chains = "TRB"
    ),
    annotations = annotations
  )
}


test_that("ImmGLIPH configuration and tool use the common interface", {
  config <- conf_immgliph(
    v_gene_col = "v_call",
    patient_col = "donor",
    n_cores = 2L,
    seed = 123L
  )
  tool <- tool_immgliph(config)

  expect_s3_class(config, "ImmGliphToolConfig")
  expect_s3_class(tool, "ImmGliphTool")
  expect_equal(config$method, "gliph2")
  expect_equal(config$n_cores, 2L)
  expect_error(conf_immgliph(method = "unknown"))
  expect_error(conf_immgliph(n_cores = 0L))
})


test_that("ImmGLIPH returns lazy groups and annotates ImmunData", {
  testthat::local_mocked_bindings(
    immgliph_runtime_info = function() {
      list(installed = TRUE, version = "test", path = "/mock/immGLIPH")
    },
    immgliph_call = function(arguments) {
      expect_named(
        arguments$cdr3_sequences,
        c("CDR3b", "TRBV", "patient")
      )
      list(
        cluster_list = list(
          group_1 = data.frame(CDR3b = c("CASSAAA", "CASSAAB"))
        ),
        connections = data.frame(
          source = "CASSAAA",
          target = "CASSAAB"
        ),
        cluster_properties = data.frame(
          group_id = "group_1",
          group_size = 2L
        )
      )
    },
    .package = "immundapter"
  )

  idata <- make_immgliph_test_idata()
  tool <- tool_immgliph(conf_immgliph(
    v_gene_col = "v_call",
    patient_col = "donor",
    seed = 123L
  ))
  run <- tool$run(idata, output_dir = withr::local_tempdir())

  expect_s3_class(run, "ImmGliphRun")
  expect_equal(run$tool, "immgliph")
  expect_setequal(
    names(run$artifacts),
    c(
      "result",
      "input",
      "mapping",
      "groups",
      "connections",
      "cluster_properties"
    )
  )
  expect_true(all(vapply(
    run$artifacts,
    fs::file_exists,
    logical(1L)
  )))

  results <- run$results()
  expect_s3_class(results, "ImmGliphResults")
  expect_identical(results$run, run)
  expect_s3_class(results$groups, "duckplyr_df")
  expect_s3_class(results$connections, "duckplyr_df")
  expect_s3_class(results$cluster_properties, "duckplyr_df")
  expect_s3_class(results$mapping, "duckplyr_df")
  expect_identical(results$raw_output, run$artifacts$result)

  groups <- results$groups |>
    dplyr::arrange(imd_receptor_id) |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(groups), 2L)
  expect_equal(
    groups$immgliph_group_id,
    c("group_1", "group_1")
  )
  expect_equal(groups$junction_aa, c("CASSAAA", "CASSAAB"))
  expect_equal(nrow(dplyr::collect(results$connections)), 1L)
  expect_equal(nrow(dplyr::collect(results$cluster_properties)), 1L)

  annotated <- results$annotate(idata)
  annotated_groups <- annotated$annotations |>
    dplyr::arrange(imd_receptor_id) |>
    dplyr::select(imd_receptor_id, immgliph_group_id) |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(
    annotated_groups$immgliph_group_id,
    c("group_1", "group_1", NA_character_)
  )
})


test_that("ImmGLIPH requires ImmunData with TRB CDR3 sequences", {
  input <- withr::local_tempfile(fileext = ".tsv")
  writeLines("junction_aa", input)

  expect_error(
    tool_immgliph()$run(input, output_dir = withr::local_tempdir()),
    "ImmunData",
    fixed = TRUE
  )

  idata <- make_immgliph_test_idata()
  expect_error(
    immgliph_prepare_immundata(
      idata,
      conf_immgliph(cdr3_col = "missing_cdr3")
    ),
    "missing required column(s): missing_cdr3",
    fixed = TRUE
  )
})
