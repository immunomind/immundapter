make_test_tool <- function(name = "test-tool") {
  TestToolConfig <- R6::R6Class(
    classname = paste0("TestToolConfig", sample.int(1000000L, 1L)),
    inherit = AbstractConfig,
    cloneable = FALSE,
    public = list(
      value = NULL,
      initialize = function(value = "default") {
        checkmate::assert_string(value)
        super$initialize(list(value = value))
      }
    )
  )

  config_class <- TestToolConfig
  TestTool <- R6::R6Class(
    classname = paste0("TestTool", sample.int(1000000L, 1L)),
    inherit = AbstractTool,
    cloneable = FALSE,
    public = list(
      initialize = function(config = config_class$new()) {
        super$initialize(name, config)
      }
    ),
    private = list(
      help_impl = function() {
        "test help"
      },
      validate_input_impl = function(input, config) {
        if (inherits(input, "ImmunData")) {
          return(input)
        }
        checkmate::assert_string(input, min.chars = 1L)
        checkmate::assert_file_exists(input)
        fs::path_abs(input)
      },
      run_impl = function(input, config, run_dir) {
        list(run_dir = run_dir)
      }
    )
  )

  list(
    tool = TestTool$new(),
    config_class = TestToolConfig
  )
}

test_that("AbstractTool exposes help instead of install", {
  tool <- make_test_tool()$tool

  expect_identical(tool$help(), "test help")
  expect_true("help" %in% names(tool))
  expect_false("install" %in% names(tool))
})

test_that("AbstractTool runs with a file and explicit output directory", {
  root <- withr::local_tempdir()
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  output_dir <- fs::path(root, "results")
  fixture <- make_test_tool()

  run <- fixture$tool$run(input, output_dir = output_dir)

  expect_s3_class(run, "ToolRun")
  expect_equal(run$status, "success")
  expect_equal(run$input$type, "files")
  expect_equal(unname(run$input$paths), as.character(fs::path_abs(input)))
  expect_equal(
    as.character(fs::path_dir(run$path)),
    as.character(fs::path_abs(output_dir))
  )
  expect_true(fs::dir_exists(run$path))
})

test_that("AbstractTool defaults file runs directly to the user home", {
  root <- withr::local_tempdir()
  home <- fs::path(root, "home")
  fs::dir_create(home)
  withr::local_envvar(HOME = home)
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  fixture <- make_test_tool()

  run <- fixture$tool$run(input)

  expect_equal(
    as.character(fs::path_dir(run$path)),
    as.character(fs::path_abs(home))
  )
})

test_that("AbstractTool uses the ImmunData snapshot artifacts path", {
  root <- withr::local_tempdir()
  artifacts_path <- fs::path(root, "artifacts", "analysis", "v001")

  TestImmunData <- R6::R6Class(
    classname = "ImmunData",
    public = list(
      provenance = NULL,
      initialize = function(provenance) {
        self$provenance <- provenance
      }
    )
  )
  input <- TestImmunData$new(list(
    current_path = fs::path(root, "snapshots", "analysis", "v001"),
    artifacts_path = artifacts_path,
    snapshot_id = "snapshot-1"
  ))
  fixture <- make_test_tool()

  run <- fixture$tool$run(input)

  expect_equal(
    as.character(fs::path_dir(run$path)),
    as.character(fs::path_abs(artifacts_path))
  )
  expect_equal(run$input$type, "immundata")
  expect_equal(run$input$snapshot_id, "snapshot-1")
})

test_that("AbstractTool requires output_dir for unpersisted ImmunData", {
  TestImmunData <- R6::R6Class(
    classname = "ImmunData",
    public = list(
      provenance = list(
        current_path = NULL,
        artifacts_path = NULL,
        snapshot_id = NULL
      )
    )
  )
  fixture <- make_test_tool()

  expect_error(
    fixture$tool$run(TestImmunData$new()),
    "without a persisted snapshot artifacts path",
    fixed = TRUE
  )
})

test_that("AbstractTool uses an immutable per-run configuration", {
  root <- withr::local_tempdir()
  input <- fs::path(root, "input.tsv")
  writeLines("value", input)
  fixture <- make_test_tool()
  config <- fixture$config_class$new("custom")

  run <- fixture$tool$run(
    input,
    config = config,
    output_dir = fs::path(root, "results")
  )

  expect_equal(run$config$value, "custom")
  expect_error(
    config$value <- "changed-after-run",
    "locked binding",
    fixed = TRUE
  )
})
