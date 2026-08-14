test_that("IgBLAST runs ImmunData through preserved AIRR artifacts", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = function() {
      list(
        package_installed = TRUE,
        igblastr_version = "1.0.23",
        igblast_installed = TRUE,
        igblast_version = "1.22.0",
        igblast_build = "igblast 1.22.0, mock build",
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
    c("input", "config", "mapping", "query_status", "airr")
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
  expect_equal(manifest$runtime$igblastr_version, "1.0.23")
  expect_equal(manifest$runtime$igblast_version, "1.22.0")
  expect_equal(
    manifest$runtime$igblast_build,
    "igblast 1.22.0, mock build"
  )
  expect_equal(manifest$runtime$igblast_root, "/mock/igblast")
  expect_equal(manifest$number_original_inputs, 3L)
  expect_equal(manifest$number_submitted_inputs, 2L)
  expect_equal(manifest$number_unique_submitted_sequences, 1L)
  expect_equal(manifest$result_status_counts$assigned, 2L)
  expect_equal(manifest$result_status_counts$unassigned, 0L)
  expect_equal(manifest$result_status_counts$missing_output, 0L)
  expect_equal(manifest$result_status_counts$not_submitted, 1L)
  expect_equal(manifest$effective_args$num_threads, 3L)
  expect_equal(manifest$effective_args$num_alignments_V, 1L)
  expect_false(any(c("query", "out") %in% names(manifest$effective_args)))

  results <- run$results()
  expect_s3_class(results, "IgblastResults")
  expect_identical(results$run, run)
  expect_s3_class(results$airr, "duckplyr_df")
  expect_s3_class(results$mapping, "duckplyr_df")
  expect_s3_class(results$query_status, "duckplyr_df")
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

  query_status <- results$query_status |>
    dplyr::collect() |>
    as.data.frame()
  expect_identical(
    colnames(query_status),
    c(
      "input_id",
      "input_id_type",
      "igblast_query_id",
      "submission_status",
      "result_status",
      "reason"
    )
  )
  expect_equal(nrow(query_status), 3L)
  expect_true(all(
    query_status$result_status[
      query_status$input_id %in% c("chain-1", "chain-2")
    ] == "assigned"
  ))
  expect_equal(
    query_status$result_status[query_status$input_id == "chain-3"],
    "not_submitted"
  )
  expect_equal(
    query_status$reason[query_status$input_id == "chain-3"],
    "missing_sequence"
  )
  expect_error(
    results$query_status <- NULL,
    "`query_status` is read-only",
    fixed = TRUE
  )

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

test_that("IgBLAST runs FASTA input through preserved AIRR artifacts", {
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
      query_ids <- sub("^>", "", fasta[seq.int(
        1L,
        length(fasta),
        by = 2L
      )])
      sequences <- fasta[seq.int(2L, length(fasta), by = 2L)]

      utils::write.table(
        data.frame(
          sequence_id = query_ids,
          sequence = sequences,
          sequence_alignment = sequences,
          germline_alignment = sequences
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

  input <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(
    ">query-2 full description",
    "acgt",
    "nnnn",
    ">query-1",
    "ACGTNNNN",
    ">query-3",
    "ttgg"
  ), input)

  run <- tool_igblast()$run(
    input,
    output_dir = withr::local_tempdir()
  )

  expect_equal(run$input$type, "files")
  expect_equal(unname(run$input$paths), as.character(fs::path_abs(input)))
  expect_equal(
    readLines(run$artifacts$input),
    c(
      ">query-2",
      "ACGTNNNN",
      ">query-1",
      "ACGTNNNN",
      ">query-3",
      "TTGG"
    )
  )

  manifest <- jsonlite::read_json(
    run$artifacts$config,
    simplifyVector = TRUE
  )
  expect_equal(manifest$input_type, "fasta")

  results <- run$results()
  mapping <- results$mapping |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(mapping), 3L)
  expect_setequal(mapping$fasta_id, c("query-1", "query-2", "query-3"))
  expect_equal(
    mapping$fasta_header[mapping$fasta_id == "query-2"],
    "query-2 full description"
  )
  expect_equal(
    mapping$igblast_query_id,
    mapping$fasta_id
  )

  airr <- results$airr |>
    dplyr::collect() |>
    as.data.frame()
  expect_equal(nrow(airr), 3L)
  expect_equal(airr$sequence_id, c("query-2", "query-1", "query-3"))

  expect_error(
    results$annotate(make_igblast_test_immundata()),
    "results produced from FASTA input cannot annotate",
    fixed = TRUE
  )
  expect_error(
    results$final,
    "No final merged table is available",
    fixed = TRUE
  )
})

test_that("IgBLAST preserves numeric-looking identifiers as text", {
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
      query_ids <- sub("^>", "", fasta[seq.int(
        1L,
        length(fasta),
        by = 2L
      )])
      writeLines(
        c(
          "sequence_id\tv_call\tj_call",
          paste(query_ids, "TRBV1*01", "TRBJ1*01", sep = "\t")
        ),
        arguments$out
      )
      invisible(NULL)
    },
    .package = "immundapter"
  )

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(">001", "ACGT", ">2", "TGCA"), fasta)
  annotations <- withr::local_tempfile(fileext = ".csv")
  writeLines(c(
    "contig_id,barcode,reads",
    "001,cell-1,100",
    "2,cell-2,80"
  ), annotations)

  results <- tool_igblast()$run(
    c(fasta = fasta, annotations = annotations),
    output_dir = withr::local_tempdir()
  )$results()
  airr <- results$airr |>
    dplyr::collect() |>
    as.data.frame()
  final <- results$final |>
    dplyr::collect() |>
    as.data.frame()

  expect_type(airr$sequence_id, "character")
  expect_equal(airr$sequence_id, c("001", "2"))
  expect_type(final$sequence_id, "character")
  expect_equal(final$sequence_id, c("001", "2"))
  expect_equal(final$contig_id, c("001", "2"))
  expect_type(final$reads, "double")
})

test_that("IgBLAST merges FASTA results with 10x annotations", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = function() {
      list(
        package_installed = TRUE,
        igblastr_version = "1.0.23",
        igblast_installed = TRUE,
        igblast_version = "1.22.0",
        igblast_root = "/mock/igblast",
        germline_db = "IMGT.mock.Homo_sapiens.IGH+IGK+IGL",
        c_region_db = "IMGT.mock.Homo_sapiens.C"
      )
    },
    igblast_call = function(arguments) {
      fasta <- readLines(arguments$query)
      query_ids <- sub("^>", "", fasta[seq.int(
        1L,
        length(fasta),
        by = 2L
      )])
      sequences <- fasta[seq.int(2L, length(fasta), by = 2L)]

      utils::write.table(
        data.frame(
          sequence_id = query_ids,
          sequence = sequences,
          locus = c("IGH", "IGK"),
          productive = TRUE,
          v_call = c("IGHV1-2*01", "IGKV1-5*01"),
          d_call = c("IGHD1-1*01", NA_character_),
          j_call = c("IGHJ4*01", "IGKJ1*01"),
          c_call = c("IGHM*01", "IGKC*01"),
          cdr3 = c("TGTGCC", "TGTAAA"),
          junction = c("TGTGCCAGC", "TGTAAATTT"),
          junction_aa = c("CAS", "CKF")
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

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(
    ">cell-1_contig-1",
    "ACGT",
    ">cell-1_contig-2",
    "TGCA"
  ), fasta)
  annotations <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      contig_id = c("cell-1_contig-1", "cell-1_contig-2"),
      barcode = c("cell-1", "cell-1"),
      chain = c("IGH", "IGK"),
      v_gene = c("OLD_IGHV", "OLD_IGKV"),
      d_gene = c("OLD_IGHD", ""),
      j_gene = c("OLD_IGHJ", "OLD_IGKJ"),
      c_gene = c("OLD_IGHM", "OLD_IGKC"),
      productive = FALSE,
      cdr3 = c("OLD_AA_H", "OLD_AA_K"),
      cdr3_nt = c("OLD_NT_H", "OLD_NT_K"),
      reads = c(100L, 80L),
      umis = c(10L, 8L),
      custom_metadata = c("heavy", "light")
    ),
    annotations,
    row.names = FALSE
  )

  run <- tool_igblast()$run(
    c(fasta = fasta, annotations = annotations),
    output_dir = withr::local_tempdir()
  )

  expect_equal(run$input$type, "files")
  expect_setequal(
    names(run$artifacts),
    c("input", "config", "mapping", "query_status", "airr", "final")
  )
  manifest <- jsonlite::read_json(
    run$artifacts$config,
    simplifyVector = TRUE
  )
  expect_equal(manifest$input_type, "fasta_annotations")
  expect_equal(
    manifest$runtime$c_region_db,
    "IMGT.mock.Homo_sapiens.C"
  )

  results <- run$results()
  expect_s3_class(results$final, "duckplyr_df")
  expect_identical(results$final_output, run$artifacts$final)
  final <- results$final |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(final$sequence_id, final$contig_id)
  expect_equal(final$barcode, c("cell-1", "cell-1"))
  expect_equal(final$reads, c(100L, 80L))
  expect_equal(final$umis, c(10L, 8L))
  expect_equal(final$custom_metadata, c("heavy", "light"))
  expect_equal(final$locus, c("IGH", "IGK"))
  expect_equal(final$v_call, c("IGHV1-2*01", "IGKV1-5*01"))
  expect_equal(final$result_status, c("assigned", "assigned"))
  expect_true(all(final$productive))
  expect_equal(final$old_chain, c("IGH", "IGK"))
  expect_equal(final$old_v_gene, c("OLD_IGHV", "OLD_IGKV"))
  expect_equal(final$old_c_gene, c("OLD_IGHM", "OLD_IGKC"))
  expect_equal(final$old_productive, c(FALSE, FALSE))
  expect_equal(final$old_cdr3, c("OLD_AA_H", "OLD_AA_K"))
  expect_equal(final$old_cdr3_nt, c("OLD_NT_H", "OLD_NT_K"))
  expect_false(any(c(
    "chain", "v_gene", "d_gene", "j_gene", "c_gene", "cdr3_nt"
  ) %in% colnames(final)))

  raw_airr <- results$airr |>
    dplyr::collect() |>
    as.data.frame()
  expect_false("barcode" %in% colnames(raw_airr))

  idata <- immundata::read_repertoires(
    path = results$final_output,
    schema = immundata::make_receptor_schema(
      features = c("v_call", "j_call", "junction_aa"),
      chains = c("IGH", "IGK|IGL")
    ),
    barcode_col = "barcode",
    locus_col = "locus",
    umi_col = "umis",
    preprocess = immundata::make_default_preprocessing("10x"),
    postprocess = NULL,
    rename_columns = NULL,
    output_folder = withr::local_tempdir()
  )
  clean_annotations <- idata$annotations |>
    dplyr::collect() |>
    as.data.frame()
  clean_receptors <- idata$receptors |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(clean_annotations), 2L)
  expect_equal(nrow(clean_receptors), 1L)
  expect_setequal(clean_annotations$locus, c("IGH", "IGK"))
  expect_equal(clean_annotations$imd_barcode, c("cell-1", "cell-1"))
  expect_true(all(c(
    "sequence_id", "v_call", "j_call", "junction_aa", "old_v_gene"
  ) %in% colnames(clean_annotations)))
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

  fixture <- read_igblast_fixture()$annotations |>
    dplyr::select(source_id, sequence) |>
    dplyr::collect() |>
    as.data.frame()
  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(
    as.vector(rbind(paste0(">", fixture$source_id), fixture$sequence)),
    fasta
  )

  fasta_run <- tool_igblast()$run(
    fasta,
    output_dir = withr::local_tempdir()
  )
  fasta_airr <- fasta_run$results()$airr |>
    dplyr::select(sequence_id, locus, v_call, j_call, productive) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(fasta_airr), 3L)
  expect_setequal(fasta_airr$sequence_id, fixture$source_id)
  expect_true(all(fasta_airr$locus == "TRB"))
  expect_true(all(grepl("^TRBV", fasta_airr$v_call)))
  expect_true(all(grepl("^TRBJ", fasta_airr$j_call)))

  annotations <- withr::local_tempfile(fileext = ".csv")
  utils::write.csv(
    data.frame(
      sequence_id = fixture$source_id,
      barcode = paste0("cell-", seq_len(nrow(fixture))),
      chain = "OLD_TRB",
      v_gene = "OLD_TRBV",
      productive = FALSE
    ),
    annotations,
    row.names = FALSE
  )
  paired_run <- tool_igblast()$run(
    c(fasta = fasta, annotations = annotations),
    output_dir = withr::local_tempdir()
  )
  paired_final <- paired_run$results()$final |>
    dplyr::collect() |>
    as.data.frame()

  expect_setequal(paired_final$sequence_id, fixture$source_id)
  expect_equal(paired_final$old_chain, rep("OLD_TRB", nrow(fixture)))
  expect_equal(paired_final$old_v_gene, rep("OLD_TRBV", nrow(fixture)))
  expect_equal(
    paired_final$productive,
    fasta_airr$productive[match(
      paired_final$sequence_id,
      fasta_airr$sequence_id
    )]
  )
})
