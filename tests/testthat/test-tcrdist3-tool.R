test_that("tcrdist3 tool clusters file and ImmunData inputs", {
  testthat::skip_on_cran()

  config <- conf_tcrdist3(
    compute_distances = FALSE,
    store_all_cdr = FALSE,
    radius = 24L
  )
  python <- Sys.which(config$python)
  if (!nzchar(python)) {
    stop(
      "The managed tcrdist3 Python environment is required outside CRAN: ",
      config$python
    )
  }
  import_status <- suppressWarnings(system2(
    python,
    c("-c", shQuote("import tcrdist")),
    stdout = FALSE,
    stderr = FALSE
  ))
  if (!identical(import_status, 0L)) {
    stop(
      "tcrdist3 must be installed in the managed environment outside CRAN."
    )
  }

  fixture <- system.file(
    "extdata",
    "tcrdist3",
    "tcr-cluster-simple.tsv",
    package = "immundapter"
  )
  snapshot <- system.file(
    "extdata",
    "tcrdist3",
    "immundata-tcr-cluster-simple",
    package = "immundapter"
  )
  paired_fixture <- system.file(
    "extdata",
    "tcrdist3",
    "tcr-cluster-paired-simple.tsv",
    package = "immundapter"
  )
  expect_true(nzchar(fixture))
  expect_true(nzchar(snapshot))
  expect_true(nzchar(paired_fixture))

  tool <- tool_tcrdist3(config)
  output_dir <- withr::local_tempdir()

  file_run <- tool$run(fixture, output_dir = output_dir)
  expect_true(inherits(file_run, "Tcrdist3Run"))
  expect_true(all(vapply(
    file_run$artifacts,
    fs::file_exists,
    logical(1L)
  )))
  file_results <- file_run$results()
  expect_true(inherits(file_results, "Tcrdist3Results"))
  expect_identical(file_results$run, file_run)
  expect_identical(file_results$distances, file_run$artifacts$distances)
  expect_error(
    file_results$distances <- "replacement.npy",
    "`distances` is read-only.",
    fixed = TRUE
  )
  file_clusters <- file_results$clusters |>
    dplyr::select(
      input_id,
      tcrdist_cluster_id,
      tcrdist_cluster_size
    ) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(file_clusters), 7L)
  expect_equal(
    sort(unique(file_clusters$tcrdist_cluster_size)),
    c(1L, 2L, 4L)
  )
  expect_equal(
    file_clusters$tcrdist_cluster_id[file_clusters$input_id == "tcr1"],
    file_clusters$tcrdist_cluster_id[file_clusters$input_id == "tcr3"]
  )
  expect_false(
    file_clusters$tcrdist_cluster_id[file_clusters$input_id == "tcr1"] ==
      file_clusters$tcrdist_cluster_id[file_clusters$input_id == "tcr4"]
  )
  file_edges <- file_results$edges |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(file_edges), 7L)

  file_mapping <- file_results$mapping |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(
    file_clusters$tcrdist_cluster_id[
      file_clusters$input_id == "tcr1"
    ],
    file_clusters$tcrdist_cluster_id[
      file_clusters$input_id == "tcr1-copy"
    ]
  )
  expect_equal(
    file_clusters$tcrdist_cluster_size[
      file_clusters$input_id == "tcr1-copy"
    ],
    4L
  )
  expect_equal(
    file_mapping$deduplicated_receptor_id[
      file_mapping$input_id == "tcr1"
    ],
    file_mapping$deduplicated_receptor_id[
      file_mapping$input_id == "tcr1-copy"
    ]
  )
  expect_equal(
    length(unique(file_mapping$deduplicated_receptor_id)),
    6L
  )

  idata <- immundata::read_immundata(snapshot, verbose = FALSE)
  idata_run <- tool$run(idata, output_dir = output_dir)
  idata_results <- idata_run$results()
  annotated <- idata_results$annotate(idata)
  idata_clusters <- annotated$annotations |>
    dplyr::select(
      expected_cluster,
      tcrdist_cluster_id,
      tcrdist_cluster_size
    ) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(idata_clusters), 7L)
  expect_equal(
    sort(unique(idata_clusters$tcrdist_cluster_size)),
    c(1L, 2L, 3L)
  )
  expect_equal(
    length(unique(idata_clusters$tcrdist_cluster_id)),
    length(unique(idata_clusters$expected_cluster))
  )
  expect_true(all(vapply(
    split(
      idata_clusters$tcrdist_cluster_id,
      idata_clusters$expected_cluster
    ),
    function(cluster) length(unique(cluster)) == 1L,
    logical(1L)
  )))

  paired_idata <- immundata::read_repertoires(
    path = paired_fixture,
    schema = immundata::make_receptor_schema(
      features = c("v_call", "j_call", "junction_aa"),
      chains = c("TRA", "TRB")
    ),
    barcode_col = "cell_id",
    locus_col = "locus",
    umi_col = "duplicate_count",
    output_folder = fs::path(output_dir, "paired-immundata"),
    repertoire_schema = NULL,
    preprocess = NULL,
    postprocess = NULL,
    rename_columns = NULL
  )
  paired_config <- conf_tcrdist3(
    chains = c("alpha", "beta"),
    cdr3_cols = c(alpha = "junction_aa", beta = "junction_aa"),
    v_gene_cols = c(alpha = "v_call", beta = "v_call"),
    j_gene_cols = c(alpha = "j_call", beta = "j_call"),
    compute_distances = FALSE,
    store_all_cdr = FALSE,
    radius = 24L
  )
  paired_tool <- tool_tcrdist3(paired_config)
  paired_run <- paired_tool$run(paired_idata, output_dir = output_dir)
  paired_results <- paired_run$results()
  paired_mapping <- paired_results$mapping |>
    dplyr::collect() |>
    as.data.frame()
  paired_annotated <- paired_results$annotate(paired_idata)
  paired_annotations_result <- paired_annotated$annotations |>
    dplyr::select(cell_id, imd_receptor_id, tcrdist_cluster_id) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(dplyr::collect(paired_idata$receptors)), 2L)
  expect_equal(nrow(dplyr::collect(paired_results$clusters)), 2L)
  expect_equal(nrow(paired_mapping), 2L)
  expect_equal(
    length(unique(paired_mapping$deduplicated_receptor_id)),
    2L
  )
  expect_equal(
    nrow(paired_annotations_result),
    6L
  )
  expect_equal(
    unique(paired_annotations_result$imd_receptor_id[
      paired_annotations_result$cell_id == "paired1"
    ]),
    unique(paired_annotations_result$imd_receptor_id[
      paired_annotations_result$cell_id == "paired1-copy"
    ])
  )
  expect_true(all(vapply(
    split(
      paired_annotations_result$tcrdist_cluster_id,
      paired_annotations_result$imd_receptor_id
    ),
    function(cluster) length(unique(cluster)) == 1L,
    logical(1L)
  )))
})
