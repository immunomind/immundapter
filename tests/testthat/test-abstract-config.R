test_that("AbstractConfig fields are immutable after validation", {
  config <- conf_tcrdist3(cpus = 2L)

  expect_equal(config$cpus, 2L)
  expect_error(
    config$cpus <- 0L,
    "locked binding",
    fixed = TRUE
  )
  expect_equal(config$cpus, 2L)
})

test_that("AbstractConfig objects cannot be cloned into mutable copies", {
  config <- conf_tcrdist3(cpus = 2L)

  expect_null(config$clone)
})
