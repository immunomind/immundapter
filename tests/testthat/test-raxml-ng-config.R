test_that("RAxML-NG configuration has family-level defaults", {
  config <- conf_raxml_ng()

  expect_equal(config$alignment_col, "clone_sequence_alignment")
  expect_equal(config$family_col, "clone_id")
  expect_equal(config$tip_col, immundata::imd_schema("receptor"))
  expect_equal(config$locus, "IGH")
  expect_equal(config$analysis, "search")
  expect_equal(config$model, "GTR+G")
  expect_null(config$seed)
  expect_null(config$threads)
  expect_false(config$ancestral)
  expect_null(config$partition_col)
  expect_null(config$brlen)
  expect_named(
    config$as_list(),
    c(
      "alignment_col", "family_col", "tip_col", "locus", "analysis",
      "model", "seed", "threads", "ancestral", "partition_col", "brlen"
    )
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
  expect_silent(conf_raxml_ng(locus = "TRA"))
  expect_error(conf_raxml_ng(alignment_col = "clone_id"))
})


test_that("RAxML-NG configuration validates reproducibility settings", {
  config <- conf_raxml_ng(model = "HKY+G", seed = 17, threads = 4)

  expect_equal(config$model, "HKY+G")
  expect_equal(config$seed, 17)
  expect_equal(config$threads, 4)
  expect_error(conf_raxml_ng(model = ""))
  expect_error(conf_raxml_ng(model = "GTR + G"), "whitespace", fixed = TRUE)
  expect_error(conf_raxml_ng(seed = 0))
  expect_error(conf_raxml_ng(seed = 1.5))
  expect_error(conf_raxml_ng(threads = 0))
  expect_error(conf_raxml_ng(threads = "auto"))
})


test_that("RAxML-NG ancestral reconstruction is opt-in", {
  expect_true(conf_raxml_ng(ancestral = TRUE)$ancestral)
  expect_error(conf_raxml_ng(ancestral = 1))
  expect_error(conf_raxml_ng(ancestral = NA))
})


test_that("RAxML-NG branch linkage is partition-only and defaults to scaled", {
  config <- conf_raxml_ng(partition_col = "partitions")

  expect_equal(config$partition_col, "partitions")
  expect_equal(config$brlen, "scaled")
  expect_equal(
    conf_raxml_ng(partition_col = "partitions", brlen = "unlinked")$brlen,
    "unlinked"
  )
  expect_error(
    conf_raxml_ng(brlen = "scaled"),
    "requires `partition_col`",
    fixed = TRUE
  )
  expect_error(conf_raxml_ng(partition_col = "partitions", brlen = "other"))
})
