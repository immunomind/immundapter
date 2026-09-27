test_that("IgBLAST runs FASTA input through preserved AIRR artifacts", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = igblast_test_runtime_info,
    igblast_call = function(arguments) {
      queries <- read_igblast_test_queries(arguments)

      write_igblast_test_airr(
        arguments,
        transform(
          queries,
          sequence_alignment = sequence,
          germline_alignment = sequence
        )
      )
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

  run <- tool_igblast(conf_igblast(
    num_threads = 3L,
    run_args = list(num_alignments_V = 1L)
  ))$run(
    input,
    output_dir = withr::local_tempdir()
  )

  expect_s3_class(run, "IgblastRun")
  expect_equal(run$tool, "igblast")
  expect_setequal(names(run$artifacts), c("input", "config", "mapping", "airr"))
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
  expect_false("sequence_col" %in% names(manifest))
  expect_equal(manifest$number_input_sequences, 3L)
  expect_equal(manifest$number_unique_sequences, 3L)
  expect_equal(manifest$fixed$outfmt, "AIRR")
  expect_equal(manifest$effective_args$num_threads, 3L)
  expect_equal(manifest$effective_args$num_alignments_V, 1L)

  results <- run$results()
  expect_s3_class(results, "IgblastResults")
  expect_s3_class(results$airr, "duckplyr_df")
  expect_s3_class(results$mapping, "duckplyr_df")
  expect_identical(results$raw_output, run$artifacts$airr)
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
    results$final,
    "No final merged table is available",
    fixed = TRUE
  )
})

test_that("IgBLAST preserves numeric-looking identifiers as text", {
  testthat::local_mocked_bindings(
    igblast_runtime_info = igblast_test_runtime_info,
    igblast_call = function(arguments) {
      queries <- read_igblast_test_queries(arguments)
      writeLines(
        c(
          "sequence_id\tv_call\tj_call",
          paste(
            queries$sequence_id,
            "TRBV1*01",
            "TRBJ1*01",
            sep = "\t"
          )
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
      igblast_test_runtime_info(
        germline_db = "IMGT.mock.Homo_sapiens.IGH+IGK+IGL",
        c_region_db = "IMGT.mock.Homo_sapiens.C"
      )
    },
    igblast_call = function(arguments) {
      queries <- read_igblast_test_queries(arguments)

      write_igblast_test_airr(
        arguments,
        data.frame(
          sequence_id = queries$sequence_id,
          sequence = queries$sequence,
          sequence_alignment = queries$sequence,
          germline_alignment = queries$sequence,
          v_germline_start = c(1L, 1L),
          v_germline_end = c(4L, 4L),
          j_germline_start = c(1L, 1L),
          j_germline_end = c(2L, 2L),
          locus = c("IGH", "IGK"),
          productive = TRUE,
          v_call = c("IGHV1-2*01", "IGKV1-5*01"),
          d_call = c("IGHD1-1*01", NA_character_),
          j_call = c("IGHJ4*01", "IGKJ1*01"),
          c_call = c("IGHM*01", "IGKC*01"),
          cdr3 = c("TGTGCC", "TGTAAA"),
          junction = c("TGTGCCAGC", "TGTAAATTT"),
          junction_aa = c("CAS", "CKF")
        )
      )
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
    c("input", "config", "mapping", "airr", "final")
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
    "sequence_id", "sequence_alignment", "v_germline_start",
    "v_germline_end", "j_germline_start", "j_germline_end",
    "v_call", "j_call", "junction_aa", "old_v_gene"
  ) %in% colnames(clean_annotations)))
  expect_equal(clean_annotations$sequence_alignment, c("ACGT", "TGCA"))
  expect_equal(clean_annotations$v_germline_end, c(4L, 4L))
})

test_that("IgBLAST runs the human TCR FASTA fixture with igblastr", {
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

  fixture <- read_igblast_fixture() |>
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

test_that("IgBLAST normalization matches Change-O 1.3.5", {
  input <- data.frame(
    fixture_case = c("ordinary", "partial_v", "lambda_insertion"),
    sequence_id = c(
      "T4_AAACCCATCCCACCTA-1_contig_2",
      "T6_ACAAGACAGGTGTATT-1_contig_1",
      "T4_CATCTTCAGACAGAAC-1_contig_2"
    ),
    locus = c("IGH", "IGK", "IGL"),
    v_call = c("IGHV3-15*01", "IGKV3-15*01", "IGLV2-8*01"),
    j_call = c("IGHJ4*02", "IGKJ2*01", "IGLJ2*01,IGLJ2A*01,IGLJ3*01"),
    sequence_alignment = c(
      paste0(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTAAAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAACGCCTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTATTAAAAGCAAAACTGATGGTGGGACAACAGACTACGCTGCACCCGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCAAAAAACACGCTGTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACCACAGATCCGAGTGGTTATTACGCAGAACTGGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
      ),
      paste0(
        "GTGATGACGCAGTCTCCAGCCACCCTGTCTGTGTCTCCAGGGGAAAGAGCCACCCTCTCTTGTAGGGCCAGTCAGAGTGTTAACACCAACTTAGCCTGGTACCAGCAAAAACCTGGCCAGGCTCCCAGGCTCCTCATCTCTGGTGCATCAATCAGGGCCACTGGTATCCCAGCCAGATTCAGTGGCAGTGGGTCTGGGACAGACTTCGCTCTCACCATCAGCAGCCTGCAGTCTGAAGATTTTGCAGTTTATTACTGTCAGCAGTATAATAAGTGGCCCCCGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC"
      ),
      paste0(
        "CAGTCTGCCCTGACTCAGCCTCCCTCCGCGTCCGGGTCTCCTGGACAGTCAGTCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGAGGTCAGTAAGCGGCGGCCCTCAGGGGTCCCTGATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCGTCTCTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCAGCTCATATGCAGGCAGCAACAATTTTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG"
      )
    ),
    germline_alignment = c(
      paste0(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGAGGCTTGGTAAAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTCAGTAACGCCTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTATTAAAAGCAAAACTGATGGTGGGACAACAGACTACGCTGCACCCGTGAAAGGCAGATTCACCATCTCAAGAGATGATTCAAAAAACACGCTGTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACCACAGANNNNAGTGGTTATTACNNNNNNNNNGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
      ),
      paste0(
        "GTGATGACGCAGTCTCCAGCCACCCTGTCTGTGTCTCCAGGGGAAAGAGCCACCCTCTCCTGCAGGGCCAGTCAGAGTGTTAGCAGCAACTTAGCCTGGTACCAGCAGAAACCTGGCCAGGCTCCCAGGCTCCTCATCTATGGTGCATCCACCAGGGCCACTGGTATCCCAGCCAGGTTCAGTGGCAGTGGGTCTGGGACAGAGTTCACTCTCACCATCAGCAGCCTGCAGTCTGAAGATTTTGCAGTTTATTACTGTCAGCAGTATAATAACTGGCCTCCGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC"
      ),
      paste0(
        "CAGTCTGCCCTGACTCAGCCTCCCTCCGCGTCCGGGTCTCCTGGACAGTCAGTCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGTGGTTATAACTATGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGAGGTCAGTAA---GCGGCCCTCAGGGGTCCCTGATCGCTTCTCTGGCTCCAAGTCTGGCAACACGGCCTCCCTGACCGTCTCTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCAGCTCATATGCAGGCAGCAACAATTTTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG"
      )
    ),
    v_germline_start = c(1L, 7L, 1L),
    v_germline_end = c(302L, 287L, 296L),
    j_germline_start = c(9L, 2L, 1L),
    j_germline_end = c(48L, 39L, 38L)
  )

  expected <- data.frame(
    fixture_case = input$fixture_case,
    sequence_id = input$sequence_id,
    sequence_alignment = c(
      paste0(
        "GAGGTGCAGCTGGTGGAGTCTGGGGGA...GGCTTGGTAAAGCCTGGGGGGTCCCTTAGACTCTCCTGTGCAGCCTCTGGATTCACTTTC............AGTAACGCCTGGATGAGCTGGGTCCGCCAGGCTCCAGGGAAGGGGCTGGAGTGGGTTGGCCGTATTAAAAGCAAAACTGATGGTGGGACAACAGACTACGCTGCACCCGTGAAA...GGCAGATTCACCATCTCAAGAGATGATTCAAAAAACACGCTGTATCTGCAAATGAACAGCCTGAAAACCGAGGACACAGCCGTGTATTACTGTACCACAGATCCGAGTGGTTATTACGCAGAACTGGACTACTGGGGCCAGGGAACCCTGGTCACCGTCTCCTCAG"
      ),
      paste0(
        "......GTGATGACGCAGTCTCCAGCCACCCTGTCTGTGTCTCCAGGGGAAAGAGCCACCCTCTCTTGTAGGGCCAGTCAGAGTGTT..................AACACCAACTTAGCCTGGTACCAGCAAAAACCTGGCCAGGCTCCCAGGCTCCTCATCTCTGGTGCA.....................TCAATCAGGGCCACTGGTATCCCA...GCCAGATTCAGTGGCAGTGGG......TCTGGGACAGACTTCGCTCTCACCATCAGCAGCCTGCAGTCTGAAGATTTTGCAGTTTATTACTGTCAGCAGTATAATAAGTGGCCCCCGTACACTTTTGGCCAGGGGACCAAGCTGGAGATCAAAC"
      ),
      paste0(
        "CAGTCTGCCCTGACTCAGCCTCCCTCC...GCGTCCGGGTCTCCTGGACAGTCAGTCACCATCTCCTGCACTGGAACCAGCAGTGACGTTGGT.........GGTTATAACTATGTCTCCTGGTACCAACAGCACCCAGGCAAAGCCCCCAAACTCATGATTTATGAGGTC.....................AGTAAGCGGCCCTCAGGGGTCCCT...GATCGCTTCTCTGGCTCCAAG......TCTGGCAACACGGCCTCCCTGACCGTCTCTGGGCTCCAGGCTGAGGATGAGGCTGATTATTACTGCAGCTCATATGCAGGCAGCAACAATTTTGTGGTATTCGGCGGAGGGACCAAGCTGACCGTCCTAG"
      )
    ),
    v_germline_start = c(1L, 1L, 1L),
    v_germline_end = c(320L, 335L, 338L),
    j_germline_start = input$j_germline_start,
    j_germline_end = input$j_germline_end
  )

  fasta <- withr::local_tempfile(fileext = ".fasta")
  writeLines(c(
    ">X92216|IGHV3-15*01|Homo sapiens|F|V-REGION|162..463|302 nt|1| | | | |302+18=320| | |",
    "gaggtgcagctggtggagtctggggga...ggcttggtaaagcctggggggtcccttagactctcctgtgcagcctctgg",
    "attcactttc............agtaacgcctggatgagctgggtccgccaggctccagggaaggggctggagtgggttg",
    "gccgtattaaaagcaaaactgatggtgggacaacagactacgctgcacccgtgaaa...ggcagattcaccatctcaaga",
    "gatgattcaaaaaacacgctgtatctgcaaatgaacagcctgaaaaccgaggacacagccgtgtattactgtaccacaga",
    ">IMGT000312|IGLV2-8*01|Homo sapiens|F|V-REGION|813705..814001|297 nt|1| | || |297+42=339| | |",
    "cagtctgccctgactcagcctccctcc...gcgtccgggtctcctggacagtcagtcaccatctcctgcactggaaccag",
    "cagtgacgttggt.........ggttataactatgtctcctggtaccaacagcacccaggcaaagcccccaaactcatga",
    "tttatgaggtc.....................agtaagcggccctcaggggtccct...gatcgcttctctggctccaag",
    "......tctggcaacacggcctccctgaccgtctctgggctccaggctgaggatgaggctgattattactgcagctcata",
    "tgcaggcagcaacaatttc",
    ">M23090|IGKV3-15*01|Homo sapiens|F|V-REGION|905..1191|287 nt|1| | | | |287+48=335| | |",
    "gaaatagtgatgacgcagtctccagccaccctgtctgtgtctccaggggaaagagccaccctctcctgcagggccagtca",
    "gagtgtt..................agcagcaacttagcctggtaccagcagaaacctggccaggctcccaggctcctca",
    "tctatggtgca.....................tccaccagggccactggtatccca...gccaggttcagtggcagtggg",
    "......tctgggacagagttcactctcaccatcagcagcctgcagtctgaagattttgcagtttattactgtcagcagta",
    "taataactggcctcc"
  ), fasta)

  query <- withr::local_tempfile(fileext = ".fasta")
  writeLines(as.vector(rbind(
    paste0(">", input$sequence_id),
    input$sequence_alignment
  )), query)
  testthat::local_mocked_bindings(
    igblast_runtime_info = igblast_test_runtime_info,
    igblast_find_imgt_v_references = function(germline_db) fasta,
    igblast_call = function(arguments) {
      write_igblast_test_airr(arguments, input)
    },
    .package = "immundapter"
  )

  run <- tool_igblast()$run(
    query,
    output_dir = withr::local_tempdir()
  )
  normalized <- run$results()$airr
  expect_s3_class(normalized, "duckplyr_df")

  columns <- c(
    "fixture_case", "sequence_id", "sequence_alignment",
    "v_germline_start", "v_germline_end",
    "j_germline_start", "j_germline_end"
  )
  actual <- normalized |>
    dplyr::select(dplyr::all_of(columns)) |>
    dplyr::arrange(fixture_case) |>
    dplyr::collect() |>
    as.data.frame()
  expected <- expected[order(expected$fixture_case), columns]
  rownames(expected) <- NULL

  expect_equal(actual, expected)
  expect_equal(nrow(actual), nrow(input))
  expect_equal(length(unique(actual$sequence_id)), nrow(actual))
})
