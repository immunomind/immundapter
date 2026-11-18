test_that("tcrdist3 respects IMMUNDAPTER_DATA", {
  root <- withr::local_tempdir()
  withr::local_envvar(c(IMMUNDAPTER_DATA = "~/tool data", R_FS_HOME = root))

  config <- conf_tcrdist3()
  executable <- if (.Platform$OS.type == "windows") {
    fs::path("Scripts", "python.exe")
  } else {
    fs::path("bin", "python")
  }

  expect_equal(
    config$python,
    fs::path(
      fs::path_abs(root),
      "tool data",
      "environments",
      "tcrdist3",
      executable
    )
  )
  expect_equal(config$mode, "cluster")
  expect_equal(config$radius, 24L)
})

test_that("tcrdist3 accepts an explicit distance-only mode", {
  config <- conf_tcrdist3(
    mode = "distance",
    python = "/opt/tcrdist3/bin/python"
  )

  expect_equal(config$mode, "distance")
})

test_that("tcrdist3 managed environment and setup commands use the user data directory", {
  root <- withr::local_tempdir()
  data_dir <- fs::path(root, "user data")
  withr::local_envvar(c(
    R_USER_DATA_DIR = data_dir,
    IMMUNDAPTER_DATA = NA
  ))
  withr::local_options(cli.width = 2000, cli.num_colors = 1)

  config <- conf_tcrdist3()
  executable <- if (.Platform$OS.type == "windows") {
    fs::path("Scripts", "python.exe")
  } else {
    fs::path("bin", "python")
  }

  expect_equal(
    config$python,
    fs::path(
      fs::path_abs(data_dir),
      "R",
      "immundapter",
      "environments",
      "tcrdist3",
      executable
    )
  )

  help <- paste(capture_messages(tool_tcrdist3(config)$help()), collapse = "\n")
  expect_match(help, "uv venv", fixed = TRUE)
  expect_match(
    help, paste0("uv pip install --python '", config$python, "'"), fixed = TRUE
  )
})

test_that("tcrdist3 preserves an explicitly configured Python", {
  config <- conf_tcrdist3(python = "/opt/tcrdist3/bin/python")

  expect_equal(config$python, "/opt/tcrdist3/bin/python")
})
