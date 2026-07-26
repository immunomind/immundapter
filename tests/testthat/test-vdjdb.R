test_that("VDJdb defaults under immundapter_data/data", {
  root <- withr::local_tempdir()
  home <- fs::path(root, "home")
  fs::dir_create(home)
  withr::local_envvar(c(
    HOME = home,
    IMMUNDAPTER_DATA = NA
  ))

  source <- data_vdjdb()

  expect_equal(
    source$output_dir,
    fs::path(fs::path_abs(home), "immundapter_data", "data")
  )
})

test_that("VDJdb respects the IMMUNDAPTER_DATA root", {
  root <- withr::local_tempdir()
  withr::local_envvar(IMMUNDAPTER_DATA = root)

  source <- data_vdjdb()

  expect_equal(
    source$output_dir,
    fs::path(fs::path_abs(root), "data")
  )
})

test_that("VDJdb accepts a custom output directory", {
  root <- withr::local_tempdir()
  output_dir <- fs::path(root, "databases")

  source <- data_vdjdb(output_dir = output_dir)

  expect_equal(source$output_dir, fs::path_abs(output_dir))
  expect_error(source$output_dir <- root, "read-only", fixed = TRUE)
})

test_that("VDJdb reports the direct database destination", {
  root <- withr::local_tempdir()
  source <- data_vdjdb(output_dir = root)
  expected <- fs::path(root, "vdjdb.slim.txt")

  expect_message(
    downloaded <- source$download(),
    regexp = "vdjdb\\.slim\\.txt"
  )
  expect_false(downloaded)
  expect_equal(
    fs::path(source$output_dir, "vdjdb.slim.txt"),
    expected
  )
})
