test_that("VDJdb defaults to the package cache", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(
    R_USER_CACHE_DIR = root,
    IMMUNDAPTER_DATA = NA
  ))

  source <- data_vdjdb()

  expect_equal(
    source$output_dir,
    fs::path(fs::path_abs(root), "R", "immundapter", "data")
  )
  expect_false(fs::dir_exists(source$output_dir))
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
  withr::local_envvar(IMMUNDAPTER_DATA = fs::path(root, "override"))
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

test_that("VDJdb downloads plain and gzip-compressed slim tables", {
  skip_if(!nzchar(Sys.which("zip")), "zip is required to build the test archive")
  root <- withr::local_tempdir()
  withr::local_dir(root)
  lines <- c("cdr3\tgene", "CASSA\tTRB")
  release_dir <- fs::dir_create("vdjdb-test")
  plain <- fs::path(release_dir, "vdjdb.slim.txt")
  compressed <- fs::path(release_dir, "vdjdb.slim.txt.gz")
  writeLines(lines, plain)
  withr::with_connection(list(file = gzfile(compressed, "wt")), {
    writeLines(lines, file)
  })

  archive <- fs::path(root, "release.zip")
  testthat::local_mocked_bindings(
    download.file = function(url, destfile, ...) {
      if (endsWith(url, "/latest-version.txt")) {
        writeLines("https://example.org/release.zip", destfile)
      } else {
        file.copy(archive, destfile)
      }
    },
    .package = "utils"
  )

  source <- data_vdjdb(output_dir = fs::path(root, "cache"))
  for (entry in c(plain, compressed)) {
    unlink(archive)
    expect_equal(utils::zip(archive, entry, flags = "-q"), 0L)
    expect_true(source$download(execute = TRUE, force = TRUE))
    expect_equal(readLines(fs::path(source$output_dir, "vdjdb.slim.txt")), lines)
    table <- source$load()
    expect_s3_class(table, "duckplyr_df")
    expect_equal(dplyr::collect(table)$cdr3, "CASSA")
  }
})
