test_that("tcrdist3 uses its managed immundapter_data environment", {
  root <- withr::local_tempdir()
  withr::local_envvar(IMMUNDAPTER_DATA = root)

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

test_that("tcrdist3 managed environment defaults under the user home", {
  root <- withr::local_tempdir()
  home <- fs::path(root, "home")
  fs::dir_create(home)
  withr::local_envvar(c(
    HOME = home,
    IMMUNDAPTER_DATA = NA
  ))

  config <- conf_tcrdist3()
  executable <- if (.Platform$OS.type == "windows") {
    fs::path("Scripts", "python.exe")
  } else {
    fs::path("bin", "python")
  }

  expect_equal(
    config$python,
    fs::path(
      fs::path_abs(home),
      "immundapter_data",
      "environments",
      "tcrdist3",
      executable
    )
  )
})

test_that("tcrdist3 preserves an explicitly configured Python", {
  config <- conf_tcrdist3(python = "/opt/tcrdist3/bin/python")

  expect_equal(config$python, "/opt/tcrdist3/bin/python")
})
