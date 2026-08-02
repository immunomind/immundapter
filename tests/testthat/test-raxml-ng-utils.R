make_raxml_test_idata <- function(annotations, chains = "IGH") {
  annotations$sequence <- gsub(
    "-",
    "",
    annotations$clone_sequence_alignment
  )
  annotations <- annotations |>
    duckplyr::as_duckdb_tibble()
  immundata::ImmunData$new(
    schema = immundata::make_receptor_schema(
      features = "sequence",
      chains = chains
    ),
    annotations = annotations
  )
}


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


test_that("RAxML-NG validates family alignments and tip identifiers", {
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


test_that("RAxML-NG builds mode-specific command arguments", {
  search <- raxml_ng_command_args(
    analysis = "search",
    alignment = "alignment.fasta",
    prefix = "result"
  )
  fast <- raxml_ng_command_args(
    analysis = "fast",
    alignment = "alignment.fasta",
    prefix = "result"
  )

  expect_true("--search" %in% search)
  expect_true("--fast" %in% fast)
  expect_identical(search[match("--model", search) + 1L], "GTR+G")
  expect_false(any(c("--tree", "--seed", "--threads") %in% search))
})
