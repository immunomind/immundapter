test_that("RAxML-NG prepares family alignments from annotations", {
  annotations <- data.frame(
    imd_receptor_id = c(1:4, 2L),
    clone_id = c(
      "family-b", "family-a", "family-a", "family-b", "family-a"
    ),
    clone_sequence_alignment = c(
      "ACGT", "ACGA", "ACGG", "ACGC", "ACGA"
    )
  )
  idata <- make_raxml_test_idata(annotations)

  prepared <- raxml_ng_prepare_immundata(
    idata,
    conf_raxml_ng()
  ) |>
    dplyr::arrange(family_id, tip_id) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(
    prepared,
    data.frame(
      family_id = c("family-a", "family-a", "family-b", "family-b"),
      tip_id = c(2L, 3L, 1L, 4L),
      aligned_sequence = c("ACGA", "ACGG", "ACGT", "ACGC")
    )
  )
  expect_equal(idata$schema_receptor$features, "sequence")
})


test_that("RAxML-NG appends backend-only tips without changing ImmunData", {
  annotations <- data.frame(
    imd_receptor_id = 1:4,
    clone_id = c("family-a", "family-a", "family-a", "family-b"),
    genotype_id = c("genotype-1", "genotype-1", "genotype-2", "genotype-3"),
    clone_sequence_alignment = c("ACGT", "ACGT", "ACGA", "TCGT"),
    partitions = c(
      "V=1-2;J=3-4", "V=1-2;J=3-4",
      "V=1-2;J=3-4", "V=1-3;J=4-4"
    )
  )
  idata <- make_raxml_test_idata(annotations)
  before <- idata$annotations |>
    dplyr::collect() |>
    as.data.frame()
  additional_tips <- raxml_test_additional_tips()

  prepared <- raxml_ng_collect_input(
    idata,
    conf_raxml_ng(
      tip_col = "genotype_id",
      partition_col = "partitions"
    ),
    additional_tips
  )
  after <- idata$annotations |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(nrow(prepared), 5L)
  expect_equal(
    prepared$tip_id,
    c(
      "family-a__germline", "genotype-1", "genotype-2",
      "family-b__germline", "genotype-3"
    )
  )
  expect_equal(
    unique(prepared$partition_spec[prepared$family_id == "family-a"]),
    "V=1-2;J=3-4"
  )
  expect_equal(before, after)
})


test_that("RAxML-NG validates additional tip families and conflicts", {
  annotations <- raxml_test_annotations()
  annotations$partitions <- "V=1-2;J=3-4"
  idata <- make_raxml_test_idata(annotations)
  config <- conf_raxml_ng(partition_col = "partitions")

  unknown_family <- data.frame(
    family_id = "family-b",
    tip_id = "germline",
    aligned_sequence = "ACGG"
  ) |>
    duckplyr::as_duckdb_tibble()
  conflicting_partition <- data.frame(
    family_id = "family-a",
    tip_id = "germline",
    aligned_sequence = "ACGG",
    partition_spec = "V=1-3;J=4-4"
  ) |>
    duckplyr::as_duckdb_tibble()
  conflicting_tip <- data.frame(
    family_id = "family-a",
    tip_id = "1",
    aligned_sequence = "ACGG"
  ) |>
    duckplyr::as_duckdb_tibble()

  expect_error(
    raxml_ng_collect_input(idata, config, unknown_family),
    "unknown family/families: family-b",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_collect_input(idata, config, conflicting_partition),
    "Conflicting partition specifications",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_collect_input(idata, config, conflicting_tip),
    "tip `1` has conflicting aligned sequences",
    fixed = TRUE
  )
})


test_that("RAxML-NG selects its configured locus from annotations", {
  annotations <- data.frame(
    imd_receptor_id = rep(1:2, each = 2L),
    locus = rep(c("IGH", "IGL"), 2L),
    clone_id = c("heavy-a", "light-a", "heavy-b", "light-b"),
    clone_sequence_alignment = c("ACGT", "TTTT", "ACGA", "TTTA")
  )
  idata <- make_raxml_test_idata(
    annotations,
    chains = c("IGH", "IGL")
  )

  heavy <- raxml_ng_prepare_immundata(
    idata,
    conf_raxml_ng(locus = "IGH")
  ) |>
    dplyr::arrange(tip_id) |>
    dplyr::collect() |>
    as.data.frame()
  light <- raxml_ng_prepare_immundata(
    idata,
    conf_raxml_ng(locus = "IGL")
  ) |>
    dplyr::arrange(tip_id) |>
    dplyr::collect() |>
    as.data.frame()

  expect_equal(heavy$family_id, c("heavy-a", "heavy-b"))
  expect_equal(heavy$aligned_sequence, c("ACGT", "ACGA"))
  expect_equal(light$family_id, c("light-a", "light-b"))
  expect_equal(light$aligned_sequence, c("TTTT", "TTTA"))
})


test_that("RAxML-NG requires family alignment annotation columns", {
  annotations <- data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
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

  expect_error(
    raxml_ng_prepare_immundata(idata, conf_raxml_ng()),
    paste0(
      "`input$annotations` is missing RAxML-NG column(s): ",
      "clone_sequence_alignment"
    ),
    fixed = TRUE
  )
})


test_that("RAxML-NG validates family alignments and tip names", {
  unaligned <- make_raxml_test_idata(data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    clone_sequence_alignment = c("ACGT", "ACGTA")
  ))
  duplicate_tips <- make_raxml_test_idata(data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    sequence_id = "duplicated-tip",
    clone_sequence_alignment = c("ACGT", "ACGA")
  ))

  expect_error(
    raxml_ng_collect_input(unaligned, conf_raxml_ng()),
    "not aligned to one length",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_collect_input(
      duplicate_tips,
      conf_raxml_ng(tip_col = "sequence_id")
    ),
    "conflicting aligned sequences",
    fixed = TRUE
  )
})


test_that("RAxML-NG builds search and fast command arguments", {
  search <- raxml_ng_command_args(
    conf_raxml_ng(model = "HKY+G", seed = 11, threads = 3),
    "/tmp/alignment file.fasta",
    "/tmp/result prefix"
  )
  fast <- raxml_ng_command_args(
    conf_raxml_ng(analysis = "fast"),
    "/tmp/alignment.fasta",
    "/tmp/result"
  )

  expect_equal(search, c(
    "--search",
    "--msa", shQuote("/tmp/alignment file.fasta"),
    "--model", shQuote("HKY+G"),
    "--prefix", shQuote("/tmp/result prefix"),
    "--seed", "11",
    "--threads", "3"
  ))
  expect_equal(fast[[1L]], "--fast")
  expect_equal(fast[match("--model", fast) + 1L], shQuote("GTR+G"))
  expect_false(any(c("--seed", "--threads") %in% fast))
})


test_that("RAxML-NG parses complete non-overlapping partitions", {
  partitions <- raxml_ng_parse_partitions("IGH=1-300;IGK=301-540", 540)

  expect_equal(
    partitions,
    data.frame(
      name = c("IGH", "IGK"),
      start = c(1, 301),
      end = c(300, 540)
    )
  )
  expect_equal(
    raxml_ng_parse_partitions("FR1=1-3;CDR1=4-6", 6)$name,
    c("FR1", "CDR1")
  )
  expect_error(
    raxml_ng_parse_partitions("A=1-3;B=3-6", 6),
    "overlap",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_parse_partitions("A=1-3;B=4-7", 6),
    "exceed",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_parse_partitions("A=1-2;B=4-6", 6),
    "cover every alignment position",
    fixed = TRUE
  )
  expect_error(
    raxml_ng_parse_partitions("A=1-3;A=4-6", 6),
    "unique",
    fixed = TRUE
  )
})


test_that("RAxML-NG writes and selects a per-family model file", {
  config <- conf_raxml_ng(
    model = "GTR+G",
    partition_col = "partitions",
    brlen = "scaled"
  )
  path <- withr::local_tempfile()
  partitions <- raxml_ng_parse_partitions("IGH=1-3;IGK=4-6", 6)
  raxml_ng_write_model_file(partitions, config$model, path)
  args <- raxml_ng_command_args(
    config,
    "/tmp/alignment.fasta",
    "/tmp/result",
    path
  )

  expect_equal(
    readLines(path),
    c("GTR+G, IGH = 1-3", "GTR+G, IGK = 4-6")
  )
  expect_equal(args[match("--model", args) + 1L], shQuote(path))
  expect_equal(args[match("--brlen", args) + 1L], "scaled")
})


test_that("RAxML-NG rejects conflicting family partitions before execution", {
  annotations <- raxml_test_annotations()
  annotations$partitions <- c("V=1-2;J=3-4", "V=1-3;J=4-4")
  idata <- make_raxml_test_idata(annotations)

  expect_error(
    raxml_ng_collect_input(
      idata,
      conf_raxml_ng(partition_col = "partitions")
    ),
    "Conflicting partition specifications",
    fixed = TRUE
  )
})


test_that("RAxML-NG permits family-specific alignment lengths and partitions", {
  annotations <- data.frame(
    imd_receptor_id = 1:4,
    clone_id = rep(c("family-a", "family-b"), each = 2L),
    clone_sequence_alignment = c("ACGT", "ACGA", "ACGTAC", "ACGTAA"),
    partitions = c(
      "V=1-2;J=3-4", "V=1-2;J=3-4",
      "V=1-3;J=4-6", "V=1-3;J=4-6"
    )
  )
  idata <- make_raxml_test_idata(annotations)

  prepared <- raxml_ng_collect_input(
    idata,
    conf_raxml_ng(partition_col = "partitions")
  )
  lengths <- tapply(
    nchar(prepared$aligned_sequence),
    prepared$family_id,
    unique
  )
  expect_equal(as.integer(lengths), c(4L, 6L))
})


test_that("RAxML-NG accepts IUPAC alignment symbols and rejects other characters", {
  accepted <- make_raxml_test_idata(data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    clone_sequence_alignment = c(
      "acgtryswkmbdhvnxo?.-",
      "ACGTRYSWKMBDHVNXO?.-"
    )
  ))
  invalid <- make_raxml_test_idata(data.frame(
    imd_receptor_id = 1:2,
    clone_id = "family-a",
    clone_sequence_alignment = c("ACGT", "ACGZ")
  ))

  prepared <- raxml_ng_collect_input(accepted, conf_raxml_ng())
  expect_equal(nrow(prepared), 2L)
  expect_true(all(prepared$aligned_sequence == "ACGTRYSWKMBDHVNXO?.-"))
  expect_error(
    raxml_ng_collect_input(invalid, conf_raxml_ng()),
    "unsupported character(s): Z",
    fixed = TRUE
  )
})
