test_that("RAxML-NG configuration has family-level defaults", {
  config <- conf_raxml_ng()

  expect_equal(config$alignment_col, "clone_sequence_alignment")
  expect_equal(config$family_col, "clone_id")
  expect_equal(config$tip_col, immundata::imd_schema("receptor"))
  expect_equal(config$locus, "IGH")
  expect_equal(config$analysis, "search")
  expect_named(
    config$as_list(),
    c("alignment_col", "family_col", "tip_col", "locus", "analysis")
  )
})

test_that("RAxML-NG family-level columns and locus are configurable", {
  config <- conf_raxml_ng(
    family_col = "family_id",
    tip_col = "sequence_id",
    locus = "IGL"
  )

  expect_equal(config$family_col, "family_id")
  expect_equal(config$tip_col, "sequence_id")
  expect_equal(config$locus, "IGL")
})

test_that("RAxML-NG configuration accepts only search and fast analyses", {
  expect_silent(conf_raxml_ng(analysis = "search"))
  expect_silent(conf_raxml_ng(analysis = "fast"))

  expect_error(conf_raxml_ng(analysis = "all"))
  expect_error(conf_raxml_ng(analysis = "bootstrap"))
  expect_error(conf_raxml_ng(analysis = "evaluate"))
  expect_error(conf_raxml_ng(analysis = "support"))
})

test_that("RAxML-NG configuration validates family-level fields", {
  expect_error(conf_raxml_ng(family_col = NULL))
  expect_error(conf_raxml_ng(tip_col = NULL))
  expect_error(conf_raxml_ng(locus = "TRA"))
  expect_error(conf_raxml_ng(alignment_col = "clone_id"))
})
